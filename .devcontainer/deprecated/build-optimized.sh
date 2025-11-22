#!/usr/bin/env bash
# build-optimized.sh - Build the multi-stage optimized Docker image
#
# This script builds the size-optimized image using multi-stage builds
# to eliminate layer bloat while preserving runtime package installation capability.
#
# Usage:
#   ./build-optimized.sh [--github-pat TOKEN] [--tag TAG]
#
# Examples:
#   ./build-optimized.sh
#   ./build-optimized.sh --github-pat ghp_xxxxx
#   ./build-optimized.sh --tag scdock-r-dev:v0.5.1

set -euo pipefail

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values (GENERIC for shareability)
GITHUB_PAT="${GITHUB_PAT:-}"
TAG="scdock-r-dev:v0.5.2"
USER_ID=1000           # Generic default (shareable image)
GROUP_ID=1000          # Generic default
USER_NAME=devuser      # Generic username
GROUP_NAME=devgroup    # Generic group
BUILD_MODE="generic"   # Track build mode for display

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --github-pat)
            GITHUB_PAT="$2"
            shift 2
            ;;
        --tag)
            TAG="$2"
            shift 2
            ;;
        --user-id)
            USER_ID="$2"
            BUILD_MODE="custom"
            shift 2
            ;;
        --group-id)
            GROUP_ID="$2"
            shift 2
            ;;
        --user)
            USER_NAME="$2"
            shift 2
            ;;
        --group)
            GROUP_NAME="$2"
            shift 2
            ;;
        --personal)
            # Convenience flag: use current user's UID/GID
            USER_ID=$(id -u)
            GROUP_ID=$(id -g)
            USER_NAME=$USER
            GROUP_NAME=$(id -gn)
            BUILD_MODE="personal"
            shift
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --github-pat TOKEN   GitHub personal access token (avoids rate limits)"
            echo "  --tag TAG            Docker image tag (default: scdock-r-dev:v0.5.1)"
            echo "  --personal           Build with YOUR UID/GID (for personal use only)"
            echo "  --user-id UID        Custom user ID (default: 1000 for shareable image)"
            echo "  --group-id GID       Custom group ID (default: 1000)"
            echo "  --user NAME          Custom username (default: devuser)"
            echo "  --group NAME         Custom group name (default: devgroup)"
            echo ""
            echo "Build Modes:"
            echo "  GENERIC (default):   Builds devuser:1000 - shareable with team/registry"
            echo "  PERSONAL:            Builds with your UID - only for your use"
            echo ""
            echo "Examples:"
            echo "  $0                              # Generic build (devuser:1000)"
            echo "  $0 --github-pat ghp_xxxxx       # Generic with PAT"
            echo "  $0 --personal                   # Personal build (your UID)"
            echo "  $0 --personal --github-pat ...  # Personal with PAT"
            exit 0
            ;;
        *)
            echo -e "${RED}Error: Unknown option $1${NC}"
            echo "Run '$0 --help' for usage"
            exit 1
            ;;
    esac
done

echo -e "${BLUE}======================================${NC}"
echo -e "${BLUE}Building Size-Optimized Docker Image${NC}"
echo -e "${BLUE}======================================${NC}"
echo ""
echo -e "${GREEN}Build Configuration:${NC}"
echo "  Mode:       $BUILD_MODE"
echo "  Tag:        $TAG"
echo "  User ID:    $USER_ID"
echo "  Group ID:   $GROUP_ID"
echo "  User:       $USER_NAME"
echo "  Group:      $GROUP_NAME"
if [ -n "$GITHUB_PAT" ]; then
    echo "  GitHub PAT: ✓ Set (first 10 chars: ${GITHUB_PAT:0:10}...)"
else
    echo -e "  GitHub PAT: ${YELLOW}✗ Not set (may hit rate limits)${NC}"
fi
echo ""
if [ "$BUILD_MODE" = "generic" ]; then
    echo -e "${GREEN}Building GENERIC image (devuser:1000) - shareable with team${NC}"
elif [ "$BUILD_MODE" = "personal" ]; then
    echo -e "${YELLOW}Building PERSONAL image (${USER_NAME}:${USER_ID}) - for your use only${NC}"
fi
echo ""

# Confirm before building
read -p "Continue with build? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Build cancelled."
    exit 0
fi

# Build start time
START_TIME=$(date +%s)

echo -e "${GREEN}Starting Docker build...${NC}"
echo ""

# Build command with BuildKit secrets (keeps GITHUB_PAT out of layers)
if [ -n "$GITHUB_PAT" ]; then
    # Write PAT to temporary file for secret mount
    PAT_FILE=$(mktemp)
    echo "$GITHUB_PAT" > "$PAT_FILE"

    DOCKER_BUILDKIT=1 docker build . \
      -f .devcontainer/Dockerfile.optimized \
      --secret id=github_pat,src="$PAT_FILE" \
      --build-arg USER_ID="$USER_ID" \
      --build-arg GROUP_ID="$GROUP_ID" \
      --build-arg USER="$USER_NAME" \
      --build-arg GROUP="$GROUP_NAME" \
      -t "$TAG" \
      2>&1 | tee build-optimized.log

    BUILD_EXIT_CODE=${PIPESTATUS[0]}

    # Clean up temp file
    rm -f "$PAT_FILE"
else
    # Build without PAT
    DOCKER_BUILDKIT=1 docker build . \
      -f .devcontainer/Dockerfile.optimized \
      --build-arg USER_ID="$USER_ID" \
      --build-arg GROUP_ID="$GROUP_ID" \
      --build-arg USER="$USER_NAME" \
      --build-arg GROUP="$GROUP_NAME" \
      -t "$TAG" \
      2>&1 | tee build-optimized.log

    BUILD_EXIT_CODE=${PIPESTATUS[0]}
fi

# Build end time
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

echo ""
if [ $BUILD_EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}======================================${NC}"
    echo -e "${GREEN}Build Completed Successfully!${NC}"
    echo -e "${GREEN}======================================${NC}"
    echo ""
    echo "  Duration: ${MINUTES}m ${SECONDS}s"
    echo "  Log:      build-optimized.log"
    echo ""

    # Show image info
    echo -e "${BLUE}Image Information:${NC}"
    docker images "$TAG" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
    echo ""

    echo -e "${YELLOW}Note: Docker may report a large 'Size' due to layer accounting.${NC}"
    echo -e "${YELLOW}Check actual filesystem usage inside the container with:${NC}"
    echo ""
    echo "  docker run --rm $TAG du -hsx /* 2>/dev/null | sort -h | tail -n 10"
    echo ""

    # Offer to extract renv.lock
    echo -e "${BLUE}Next Steps:${NC}"
    echo ""
    echo "1. Extract renv.lock (if first build):"
    echo "   CID=\$(docker create $TAG)"
    echo "   docker cp \$CID:/opt/settings/renv.lock ./renv.lock"
    echo "   docker cp \$CID:/opt/settings/R-packages-manifest.csv ./R-packages-manifest.csv"
    echo "   docker rm \$CID"
    echo ""
    echo "2. Run sanity checks:"
    echo "   docker run --rm $TAG bash -lc '.devcontainer/scripts/poststart_sanity.sh'"
    echo ""
    echo "3. Test interactively:"
    echo "   docker run --rm -it $TAG bash"
    echo ""

else
    echo -e "${RED}======================================${NC}"
    echo -e "${RED}Build Failed!${NC}"
    echo -e "${RED}======================================${NC}"
    echo ""
    echo "  Duration: ${MINUTES}m ${SECONDS}s"
    echo "  Log:      build-optimized.log"
    echo ""
    echo "Check the log file for errors."
    exit 1
fi
