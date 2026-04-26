#!/bin/bash
OBJECT_FILE=$1
FILENAME=$(basename "$OBJECT_FILE")
if [ -n "$2" ]; then
    EXPECTED_PREFIX=$2
else
    EXPECTED_PREFIX="reb_${FILENAME%.*}"
fi
OS_TYPE=$(uname -s)
EXIT_STATUS=0

echo "Expected prefix: $EXPECTED_PREFIX"
echo "------------------------------------------"

# Use nm to list symbols
FUNCTIONS=$(nm -g --defined-only "$OBJECT_FILE" | grep ' [T] ' | awk '{print $3}')

# Loop through functions and validate prefix
while read -r FUNC; do

	# Handle macOS underscore prefix
    CLEAN_FUNC="$FUNC"
    if [[ "$OS_TYPE" == "Darwin" ]]; then
        # Remove the first underscore if it exists
        CLEAN_FUNC="${FUNC#_}"
    fi
    if [[ "$CLEAN_FUNC" == "$EXPECTED_PREFIX"* ]]; then
        echo "[OK]   $CLEAN_FUNC"
    else
        echo "[FAIL] $CLEAN_FUNC (Does not match prefix '$EXPECTED_PREFIX')"
		EXIT_STATUS=1
    fi
done <<< "$FUNCTIONS"

exit $EXIT_STATUS
