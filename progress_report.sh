#!/bin/bash
# Generate a comprehensive progress report for the migration

echo "========================================================================"
echo "PYBALL MIGRATION PROGRESS REPORT"
echo "========================================================================"
echo ""
echo "Generated: $(date)"
echo ""

echo "--- OVERALL STATUS ---"
echo "Total issues:  $(bd count)"
echo "Open:          $(bd count --status open)"
echo "Closed:        $(bd count --status closed)"
echo "Ready to work: $(bd ready 2>&1 | grep -c '\[P')"
echo ""

echo "--- BY PRIORITY ---"
echo "P1 (Critical): $(bd list --priority 1 2>/dev/null | wc -l | tr -d ' ')"
echo "P2 (High):     $(bd list --priority 2 2>/dev/null | wc -l | tr -d ' ')"
echo "P3 (Normal):   $(bd list --priority 3 2>/dev/null | wc -l | tr -d ' ')"
echo ""

echo "--- PHASES STATUS ---"
for i in 1 2 3 4 5 6 7; do
    issue_id="pyball-hdg.$i"
    status=$(bd show $issue_id 2>/dev/null | grep "^Status:" | awk '{print $2}')
    title=$(bd show $issue_id 2>/dev/null | grep "^$issue_id:" | cut -d: -f2-)
    echo "Phase $i:$title - [$status]"
done
echo ""

echo "--- READY TO WORK ON ---"
bd ready | head -15
echo ""

echo "========================================================================"
echo "For detailed view: bd show pyball-hdg"
echo "For specific phase: bd show pyball-hdg.1"
echo "========================================================================"
