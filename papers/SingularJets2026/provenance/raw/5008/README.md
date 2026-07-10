# Case 5008 raw provenance

These files are copied byte-for-byte from the durable case archive at
`SBB-OhSweep-Archive-2026-07/Bo0-L14/5008`:

- `restart.sbatch`: the exact restart submission, including parent case 5004,
  snapshot `0.490625`, Slurm job name and resource request.
- `slurm-restart-sbb5008-24472158.out`: scheduler stdout.
- `slurm-restart-sbb5008-24472158.err`: scheduler stderr and the terminal
  `TIME LIMIT` record.

The restart snapshot checksum is intentionally not recorded.  The archive was
verified byte-for-byte to use case 5004 snapshot `0.490625` as case 5008's
`restart`.
