# Codex Agent Instructions

## General
- This is a Julia package intended for research and high-performance numerical computing.
- Prefer correctness, clarity, and performance over stylistic changes.
- Do not change public APIs unless explicitly instructed.

## Review Mode (default)
- Start in read-only review mode.
- First summarize package structure and intended public API.
- Identify issues grouped into: correctness, performance, API design, and maintainability.
- Do not modify files unless explicitly requested.

## Julia-Specific Guidelines
- Prioritize type stability and inference-friendly code.
- Flag allocations, abstract containers, dynamic dispatch, and non-const globals.
- Avoid over-parameterization or unnecessary use of Val/traits unless clearly beneficial.
- Assume performance-critical inner loops may exist.

## Numerical / Linear Algebra Code
- Be conservative: do not alter numerical algorithms unless a clear issue is identified.
- Flag potential numerical instability, normalization issues, or silent invariant violations.
- Prefer minimal, testable changes.

## Testing
- If suggesting changes, also suggest corresponding tests.
- Do not rewrite the test suite unless explicitly requested.

