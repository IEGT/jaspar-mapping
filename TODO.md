# TODO

## Regional motif exploration

- Provide a web service for dynamic analysis of a requested genomic interval.
  It should expose the underlying retained scores and, when requested, rescore
  the reference sequence below production storage thresholds. Parameters must
  include genome build, chromosome span, motifs, score mode, pseudocount,
  anchor threshold, local-peak radius, tandem distance, and context distance.
- Provide an accessible website that uses the service to inspect motif scores,
  TP73 regional prominence, tandem architecture, cofactor relationships,
  gene/transcript context, and available CUT&RUN evidence. Keep the deployment
  host open: `gentle.functional.domains` is one option, while a University-hosted
  site may be preferable. Results should remain downloadable in the documented
  table formats and usable by OpenClaw, Claude, and Codex through the same API.

The service is the escape hatch for focused regions: conservative genome-wide
packages need not persist every subordinate shifted alignment when the full
local score landscape can be regenerated on demand.
