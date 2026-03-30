# immunedeconv 2.1.1

## Enhancements

- `deconvolute_mcp_counter()` gains a `log_transform` parameter (`NULL`/`TRUE`/`FALSE`, default `NULL`).
  MCP-counter expects log-transformed expression data, but `deconvolute()` documents that it accepts
  raw TPM values. The new parameter resolves this mismatch:
  - `NULL` (default): auto-detects whether the input is already log-transformed (heuristic: `max(matrix) > 50`)
    and applies `log2(x + 1)` if needed, emitting an informational message.
  - `TRUE`: always applies `log2(x + 1)` transformation.
  - `FALSE`: skips transformation (use when data is already log-transformed).
  The parameter can also be passed through the generic `deconvolute()` interface via `...`.
