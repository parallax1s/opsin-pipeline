"""Real ESM2 backend for ``PLMScorer`` per spec §5 and §6.

Heavy imports (torch, transformers) happen lazily inside method calls so the
pipeline's default tests can run with neither dependency installed. The module
itself is always importable; a missing dependency only bites at ``score()``
time, with an actionable error.

Default model: ``facebook/esm2_t12_35M_UR50D`` — 35M params, ~140 MB, runs in
a few seconds on CPU for opsin-length sequences.
"""
from __future__ import annotations

from dataclasses import dataclass

from .predictions import PLMPrediction, PLMPredictionSet, sequence_sha256
from .scorer import (
    PLMScorer,
    PositionOutOfRangeError,
    STANDARD_AAS,
    UnsupportedResidueError,
    WildTypeMismatchError,
    MutationRequest,
)


class PLMBackendUnavailableError(ImportError):
    """Raised when torch / transformers isn't importable and a real PLM is requested."""


class SequenceTooLongError(ValueError):
    """Raised when the scaffold sequence exceeds the model's context window."""


class ModelUnknownError(ValueError):
    """Raised when a --model string isn't one of the supported ESM2 aliases."""


# Map friendly aliases → HuggingFace repo IDs. Extend as needed.
_MODEL_ALIASES: dict[str, str] = {
    "esm2_t6_8M_UR50D": "facebook/esm2_t6_8M_UR50D",
    "esm2_t12_35M_UR50D": "facebook/esm2_t12_35M_UR50D",
    "esm2_t30_150M_UR50D": "facebook/esm2_t30_150M_UR50D",
    "esm2_t33_650M_UR50D": "facebook/esm2_t33_650M_UR50D",
}

_DEFAULT_MAX_LENGTH = 1024


def resolve_model_id(name: str) -> str:
    if name in _MODEL_ALIASES:
        return _MODEL_ALIASES[name]
    if "/" in name:
        return name  # treat as a raw HF repo id
    raise ModelUnknownError(
        f"Unknown ESM2 model alias {name!r}. Known aliases: "
        + ", ".join(sorted(_MODEL_ALIASES))
        + ". Pass a full HuggingFace repo id (e.g. 'facebook/esm2_t6_8M_UR50D') for anything else."
    )


def _import_backend():
    """Lazy, side-effect-free attempt to import torch + transformers.

    Returns ``(torch_module, transformers_module)`` or raises
    ``PLMBackendUnavailableError`` with a useful install hint.
    """
    try:
        import torch  # noqa: F401
        import transformers  # noqa: F401
    except ImportError as exc:
        raise PLMBackendUnavailableError(
            "ESM2Scorer requires torch + transformers. Install with: "
            "pip install torch transformers. Or rerun with --mock to use "
            "the deterministic MockPLMScorer (no network, no dependencies)."
        ) from exc
    return __import__("torch"), __import__("transformers")


@dataclass
class ESM2Scorer(PLMScorer):
    """ESM2 sequence-only plausibility scorer.

    Mask-free: one forward pass on the WT sequence gives per-position log
    probabilities for every amino acid. For any (position, from_aa, to_aa),
    ``delta = log_p[pos][to_aa] - log_p[pos][from_aa]`` — O(1) per mutation
    after a single forward pass per scaffold.
    """

    model_alias: str = "esm2_t12_35M_UR50D"
    device: str = "cpu"
    max_length: int = _DEFAULT_MAX_LENGTH

    def __post_init__(self) -> None:
        # Validate the alias eagerly so `plm --model <bad>` fails fast even
        # without torch installed.
        self.hf_model_id = resolve_model_id(self.model_alias)
        self.model_id = self.model_alias
        self._tokenizer = None
        self._model = None
        self._aa_token_ids: dict[str, int] | None = None

    def score(
        self,
        *,
        scaffold_name: str,
        sequence: str,
        mutations: list[MutationRequest],
    ) -> PLMPredictionSet:
        if len(sequence) > self.max_length:
            raise SequenceTooLongError(
                f"sequence length {len(sequence)} exceeds model max_length "
                f"{self.max_length}"
            )
        self._ensure_model()
        log_probs = self._compute_log_probs(sequence)

        predictions: list[PLMPrediction] = []
        for req in mutations:
            self._validate(sequence, req)
            if req.from_aa == req.to_aa:
                delta = 0.0
            else:
                from_id = self._aa_token_ids[req.from_aa]
                to_id = self._aa_token_ids[req.to_aa]
                # log_probs shape is [seq_len, vocab]; position is 1-indexed.
                delta = float(
                    log_probs[req.position - 1, to_id] - log_probs[req.position - 1, from_id]
                )
            predictions.append(
                PLMPrediction(
                    position=req.position,
                    from_aa=req.from_aa,
                    to_aa=req.to_aa,
                    log_likelihood_delta=round(delta, 4),
                    model_id=self.model_id,
                )
            )
        return PLMPredictionSet(
            scaffold_name=scaffold_name,
            sequence_sha256=sequence_sha256(sequence),
            model_id=self.model_id,
            predictions=predictions,
        )

    # ---- internals ----

    def _ensure_model(self) -> None:
        if self._model is not None and self._tokenizer is not None:
            return
        torch, transformers = _import_backend()
        self._tokenizer = transformers.AutoTokenizer.from_pretrained(self.hf_model_id)
        self._model = transformers.AutoModelForMaskedLM.from_pretrained(self.hf_model_id)
        self._model.eval()
        if self.device != "cpu":
            self._model.to(self.device)
        # Build an amino-acid -> token-id table for fast lookup.
        aa_token_ids: dict[str, int] = {}
        for aa in STANDARD_AAS:
            tid = self._tokenizer.convert_tokens_to_ids(aa)
            if tid is None or tid == self._tokenizer.unk_token_id:
                raise RuntimeError(
                    f"ESM2 tokenizer has no token id for standard amino acid {aa!r}; "
                    "model checkpoint may be incompatible."
                )
            aa_token_ids[aa] = tid
        self._aa_token_ids = aa_token_ids

    def _compute_log_probs(self, sequence: str):
        torch, _ = _import_backend()
        encoded = self._tokenizer(
            sequence, return_tensors="pt", add_special_tokens=True
        )
        if self.device != "cpu":
            encoded = {k: v.to(self.device) for k, v in encoded.items()}
        with torch.no_grad():
            logits = self._model(**encoded).logits  # [1, L+special, vocab]
        log_probs = torch.log_softmax(logits[0], dim=-1)  # [L+special, vocab]
        # ESM2 prepends a <cls> and appends an <eos> token. Strip them so the
        # returned tensor is indexed directly by 0-based residue index.
        return log_probs[1 : 1 + len(sequence)].detach().cpu()

    def _validate(self, sequence: str, req: MutationRequest) -> None:
        if req.position < 1 or req.position > len(sequence):
            raise PositionOutOfRangeError(
                f"position {req.position} is outside sequence length {len(sequence)}"
            )
        actual_wt = sequence[req.position - 1].upper()
        if actual_wt not in STANDARD_AAS:
            raise UnsupportedResidueError(
                f"scaffold residue at position {req.position} is {actual_wt!r}; "
                "only the standard 20 amino acids are supported"
            )
        if req.to_aa not in STANDARD_AAS:
            raise UnsupportedResidueError(
                f"target residue {req.to_aa!r} at position {req.position} is not one of the standard 20"
            )
        if req.from_aa != actual_wt:
            raise WildTypeMismatchError(
                f"mutation request says from_aa={req.from_aa} at position {req.position}, "
                f"but scaffold sequence has {actual_wt}"
            )
