#!/usr/bin/env bash
# build-archr-wrapper.sh - Build the ArchR wrapper image
#
# This script builds a wrapper around the official ArchR image
# (greenleaflab/archr:1.0.3-base-r4.4) with consistent UID/GID handling
# to match the scdock-r-dev image pattern.
#
# Usage:
#   ./build-archr-wrapper.sh [OPTIONS]
#
# Examples:
#   ./build-archr-wrapper.sh                    # Generic build (devuser:1000)
#   ./build-archr-wrapper.sh --personal         # Personal build (your UID)
#   ./build-archr-wrapper.sh --tag custom:tag   # Custom tag

set -euo pipefail

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values (GENERIC for shareability)
TAG="scdock-r-archr:v0.5.1"
USER_ID=1000           # Generic default (shareable image)
GROUP_ID=1000          # Generic default
USER_NAME=devuser      # Generic username
GROUP_NAME=devgroup    # Generic group
BUILD_MODE="generic"   # Track build mode for display

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
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
            echo "  --tag TAG            Docker image tag (default: scdock-r-archr:v0.5.1)"
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
            echo "  $0 --personal                   # Personal build (your UID)"
            echo "  $0 --tag scdock-r-archr:custom  # Custom tag"
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
echo -e "${BLUE}Building ArchR Wrapper Image${NC}"
echo -e "${BLUE}======================================${NC}"
echo ""
echo -e "${GREEN}Build Configuration:${NC}"
echo "  Mode:       $BUILD_MODE"
echo "  Tag:        $TAG"
echo "  User ID:    $USER_ID"
echo "  Group ID:   $GROUP_ID"
echo "  User:       $USER_NAME"
echo "  Group:      $GROUP_NAME"
echo "  Base Image: greenleaflab/archr:1.0.3-base-r4.4"
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

# Build command
docker build . \
  -f .devcontainer/Dockerfile.archr-wrapper \
  --build-arg USER_ID="$USER_ID" \
  --build-arg GROUP_ID="$GROUP_ID" \
  --build-arg USER="$USER_NAME" \
  --build-arg GROUP="$GROUP_NAME" \
  -t "$TAG" \
  2>&1 | tee build-archr-wrapper.log

BUILD_EXIT_CODE=${PIPESTATUS[0]}

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
    echo "  Log:      build-archr-wrapper.log"
    echo ""

    # Show image info
    echo -e "${BLUE}Image Information:${NC}"
    docker images "$TAG" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
    echo ""

    echo -e "${BLUE}Next Steps:${NC}"
    echo ""
    echo "1. Test the image:"
    echo "   docker run --rm -it $TAG bash"
    echo ""
    echo "2. Use in docker-compose.yml:"
    echo "   dev-archr:"
    echo "     image: $TAG"
    echo ""
    echo "3. Update devcontainer.json:"
    echo "   \"service\": \"dev-archr\""
    echo "   \"remoteUser\": \"devuser\""
    echo "   \"updateRemoteUserUID\": true"
    echo ""

else
    echo -e "${RED}======================================${NC}"
    echo -e "${RED}Build Failed!${NC}"
    echo -e "${RED}======================================${NC}"
    echo ""
    echo "  Duration: ${MINUTES}m ${SECONDS}s"
    echo "  Log:      build-archr-wrapper.log"
    echo ""
    echo "Check the log file for errors."
    exit 1
fi
