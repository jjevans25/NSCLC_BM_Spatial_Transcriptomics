#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# marimo-check — PostToolUse hook (PROJECT_PLAN §5.6).
#
# Runs `uvx marimo check` on any marimo notebook the agent just touched. On
# failure it exits 2, which feeds stderr back to the agent so it self-corrects.
#
# A file counts as a marimo notebook only if it contains both `import marimo`
# and `@app.cell` — ordinary scripts under workflow/scripts/ are left alone.
#
# Bash is matched as well as Edit/Write: notebooks written via heredoc would
# otherwise slip past the check entirely.
# ---------------------------------------------------------------------------
set -uo pipefail

command -v uvx >/dev/null 2>&1 || exit 0

input=$(cat)
tool=$(printf '%s' "$input" | jq -r '.tool_name // empty')

candidates=$(printf '%s' "$input" | jq -r '.tool_input.file_path // empty')

if [ -z "$candidates" ] && [ "$tool" = "Bash" ]; then
  # Pull .py paths out of the command line and let the marimo-notebook test
  # below decide which of them are actually notebooks.
  candidates=$(printf '%s' "$input" \
    | jq -r '.tool_input.command // empty' \
    | grep -oE '[A-Za-z0-9_./-]+\.py' \
    | sort -u)
fi

[ -z "$candidates" ] && exit 0

status=0
report=""

while IFS= read -r file; do
  [ -z "$file" ] && continue
  case "$file" in *.py) ;; *) continue ;; esac
  [ -f "$file" ] || continue
  grep -q 'import marimo' "$file" || continue
  grep -q '@app.cell' "$file" || continue

  if ! out=$(uvx marimo check "$file" 2>&1); then
    report="${report}
--- marimo check failed: ${file} ---
${out}"
    status=2
  fi
done <<< "$candidates"

if [ "$status" -ne 0 ]; then
  printf '%s\n' "$report" >&2
  echo "" >&2
  echo "Fix only what marimo check reports, then re-save the notebook." >&2
  exit 2
fi

exit 0
