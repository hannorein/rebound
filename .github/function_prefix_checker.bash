#!/bin/bash

# Check if a file was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <file.o>"
    exit 1
fi

OBJECT_FILE=$1

# Check if the file exists
if [ ! -f "$OBJECT_FILE" ]; then
    echo "Error: File '$OBJECT_FILE' not found."
    exit 1
fi

# Get the filename without the path and extension to use as the prefix
# Example: "src/auth_utils.o" -> "auth_utils"
FILENAME=$(basename "$OBJECT_FILE")
EXPECTED_PREFIX="_reb_${FILENAME%.*}"

echo "Checking functions in: $OBJECT_FILE"
echo "Expected prefix: $EXPECTED_PREFIX"
echo "------------------------------------------"

# Use nm to list symbols
# --defined-only: ignores undefined external symbols
# grep ' [Tt] ': filters for global (T) or local (t) functions in the text section
# awk '{print $3}': extracts the symbol name (usually the 3rd column)
# -g / --extern-only: Only display symbols that are visible outside the object file
FUNCTIONS=$(nm -g --defined-only "$OBJECT_FILE" | grep ' [T] ' | awk '{print $3}')

if [ -z "$FUNCTIONS" ]; then
    echo "No functions found in $OBJECT_FILE."
    exit 0
fi

# Loop through functions and validate prefix
while read -r FUNC; do
    if [[ "$FUNC" == "$EXPECTED_PREFIX"* ]]; then
        echo "[OK]   $FUNC"
    else
        echo "[FAIL] $FUNC (Does not match prefix '$EXPECTED_PREFIX')"
    fi
done <<< "$FUNCTIONS"
