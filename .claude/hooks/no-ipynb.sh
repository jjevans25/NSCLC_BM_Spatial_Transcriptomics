#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# no-ipynb — PreToolUse guard enforcing PROJECT_PLAN §4.5:
#   "Notebooks are marimo (.py) only. .ipynb files are not created, not
#    committed, and not accepted as deliverables."
#
# PreToolUse (not PostToolUse): exit 2 here refuses the tool call, so the file
# is never written. A PostToolUse hook would only complain after the fact.
# See docs/decisions/0002-ipynb-guard-is-a-pretooluse-hook.md.
#
# Covers Write/Edit/NotebookEdit by path, and Bash by looking for commands that
# would *create* an .ipynb. `marimo convert` is explicitly allowed — that is the
# sanctioned route for Jupyter material arriving from a paper or tutorial.
# ---------------------------------------------------------------------------
set -uo pipefail

input=$(cat)

deny() {
  cat >&2 <<EOF
BLOCKED by the no-ipynb hook (PROJECT_PLAN §4.5).

  $1

This project uses marimo notebooks (.py) exclusively — they are diffable,
have no hidden state, and are runnable as scripts by a Snakemake rule.

If you are bringing in Jupyter material from a paper or tutorial, convert it:
  marimo convert <old>.ipynb -o notebooks/explore/<name>.py
Otherwise write the notebook directly as a marimo .py under notebooks/.
EOF
  exit 2
}

tool=$(printf '%s' "$input" | jq -r '.tool_name // empty')
path=$(printf '%s' "$input" | jq -r '.tool_input.file_path // .tool_input.notebook_path // empty')

case "$path" in
  *.ipynb) deny "Refused to write: $path" ;;
esac

if [ "$tool" = "Bash" ]; then
  cmd=$(printf '%s' "$input" | jq -r '.tool_input.command // empty')

  # Sanctioned conversion path — allow through.
  if printf '%s' "$cmd" | grep -Eq 'marimo[[:space:]]+convert'; then
    exit 0
  fi

  # Commands that would create or copy an .ipynb into the tree.
  if printf '%s' "$cmd" | grep -Eq '>>?[[:space:]]*"?'"'"'?[^[:space:]"'"'"']*\.ipynb'; then
    deny "Refused: shell redirection would create an .ipynb file."
  fi
  if printf '%s' "$cmd" | grep -Eq '\b(touch|cp|mv|tee|install)\b[^|;&]*\.ipynb'; then
    deny "Refused: this command would create or copy an .ipynb file."
  fi
fi

exit 0
