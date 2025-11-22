#!/usr/bin/env bash
# build.sh - Build the multi-stage optimized Docker image (restructured paths)
#
# Usage:
#   scripts/build.sh [--github-pat TOKEN] [--tag TAG] [--personal] [--user-id UID] [--group-id GID] [--user NAME] [--group NAME]
#
# Examples:
#   scripts/build.sh
#   scripts/build.sh --github-pat ghp_xxxxx
#   scripts/build.sh --tag scdock-r-dev:v0.5.2

set -euo pipefail

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

GITHUB_PAT="${GITHUB_PAT:-}"
TAG="scdock-r-dev:v0.5.2"
USER_ID=1000
GROUP_ID=1000
USER_NAME=devuser
GROUP_NAME=devgroup
BUILD_MODE="generic"

while [[ $# -gt 0 ]]; do
  case $1 in
    --github-pat) GITHUB_PAT="$2"; shift 2;;
    --tag) TAG="$2"; shift 2;;
    --user-id) USER_ID="$2"; BUILD_MODE="custom"; shift 2;;
    --group-id) GROUP_ID="$2"; shift 2;;
    --user) USER_NAME="$2"; shift 2;;
    --group) GROUP_NAME="$2"; shift 2;;
    --personal) USER_ID=$(id -u); GROUP_ID=$(id -g); USER_NAME=$USER; GROUP_NAME=$(id -gn); BUILD_MODE="personal"; shift;;
    --help)
      cat <<EOF
Usage: $0 [OPTIONS]

Options:
  --github-pat TOKEN   GitHub personal access token (avoids rate limits)
  --tag TAG            Docker image tag (default: scdock-r-dev:v0.5.2)
  --personal           Build with YOUR UID/GID (for personal use only)
  --user-id UID        Custom user ID (default: 1000)
  --group-id GID       Custom group ID (default: 1000)
  --user NAME          Custom username (default: devuser)
  --group NAME         Custom group name (default: devgroup)
EOF
      exit 0;;
    *) echo -e "${RED}Unknown option: $1${NC}"; exit 1;;
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

read -p "Continue with build? (y/N): " -n 1 -r; echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  echo "Build cancelled."; exit 0
fi

START_TIME=$(date +%s)

echo -e "${GREEN}Starting Docker build...${NC}"
echo ""

LOG_FILE=build.log

if [ -n "$GITHUB_PAT" ]; then
  PAT_FILE=$(mktemp)
  echo "$GITHUB_PAT" > "$PAT_FILE"
  DOCKER_BUILDKIT=1 docker build . \
    -f docker/base/Dockerfile \
    --secret id=github_pat,src="$PAT_FILE" \
    --build-arg USER_ID="$USER_ID" \
    --build-arg GROUP_ID="$GROUP_ID" \
    --build-arg USER="$USER_NAME" \
    --build-arg GROUP="$GROUP_NAME" \
    -t "$TAG" \
    2>&1 | tee "$LOG_FILE"
  BUILD_EXIT_CODE=${PIPESTATUS[0]}
  rm -f "$PAT_FILE"
else
  DOCKER_BUILDKIT=1 docker build . \
    -f docker/base/Dockerfile \
    --build-arg USER_ID="$USER_ID" \
    --build-arg GROUP_ID="$GROUP_ID" \
    --build-arg USER="$USER_NAME" \
    --build-arg GROUP="$GROUP_NAME" \
    -t "$TAG" \
    2>&1 | tee "$LOG_FILE"
  BUILD_EXIT_CODE=${PIPESTATUS[0]}
fi

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
  echo "  Log:      ${LOG_FILE}"
  echo ""
  echo -e "${BLUE}Image Information:${NC}"
  docker images "$TAG" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
  echo ""
  echo -e "${YELLOW}Note: Docker may report a large 'Size' due to layer accounting.${NC}"
  echo "  docker run --rm $TAG du -hsx /* 2>/dev/null | sort -h | tail -n 10"
  echo ""
  echo -e "${BLUE}Next Steps:${NC}"
  echo "1. Extract renv.lock (if first build):"
  echo "   CID=\$(docker create $TAG)"
  echo "   docker cp \$CID:/opt/settings/renv.lock ./renv.lock"
  echo "   docker cp \$CID:/opt/settings/R-packages-manifest.csv ./R-packages-manifest.csv"
  echo "   docker rm \$CID"
  echo ""
  echo "2. Run sanity checks:"
  echo "   docker run --rm $TAG bash -lc 'scripts/poststart_sanity.sh'"
  echo ""
  echo "3. Test interactively:"
  echo "   docker run --rm -it $TAG bash"
else
  echo -e "${RED}======================================${NC}"
  echo -e "${RED}Build Failed!${NC}"
  echo -e "${RED}======================================${NC}"
  echo ""
  echo "  Duration: ${MINUTES}m ${SECONDS}s"
  echo "  Log:      ${LOG_FILE}"
  echo ""
  echo "Check the log file for errors."
  exit 1
fi
