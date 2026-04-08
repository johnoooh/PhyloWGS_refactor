# Over-Splitting Audit — Go vs Python Side-by-Side

**Date:** 2026-04-08
**Branch:** go-port
**Reference Python:** /sessions/dreamy-pensive-knuth/mnt/phylowgs_cowork/phylowgs/

This document records every concrete difference between the Go port and
the upstream Python on the suspect functions identified in
docs/plans/2026-04-08-oversplitting-investigation.md.

## Format

For each function:
- Python file:line range
- Go file:line range
- Differences observed (line-by-line diff in prose)
- Severity: BLOCKER (will measurably affect over-splitting) /
            SUSPICIOUS (semantic divergence, may matter) /
            BENIGN (cosmetic / equivalent rewrite)
- Whether to fix in this plan

## Functions audited

(filled in as Phase 1 progresses)
