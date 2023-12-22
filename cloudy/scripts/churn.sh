#!/bin/sh
# Usage: call in a git checkout
#   ./churn.sh | head -20
git whatchanged | awk '/^:/ {print $6}' | sort | uniq -c | sort -nr
