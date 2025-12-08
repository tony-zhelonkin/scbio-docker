#!/usr/bin/env bash
# build-ai-enabled.sh - Build the AI-enabled scbio-docker image

set -euo pipefail

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

BASE_IMAGE_TAG="scdock-r-dev:v0.5.2"
AI_IMAGE_TAG="scdock-ai-dev:v0.5.2"

# Check if the base image exists
if ! docker image inspect "${BASE_IMAGE_TAG}" &> /dev/null; then
  echo -e "${YELLOW}Base image ${BASE_IMAGE_TAG} not found.${NC}"
  echo -e "${BLUE}Running the base image build script...${NC}"
  ./scripts/build.sh
fi

echo -e "${BLUE}======================================${NC}"
echo -e "${BLUE}Building AI-Enabled Docker Image${NC}"
echo -e "${BLUE}======================================${NC}"

DOCKER_BUILDKIT=1 docker build .
  -f docker/Dockerfile.ai-enabled
  -t "${AI_IMAGE_TAG}"

echo -e "${GREEN}======================================${NC}"
echo -e "${GREEN}AI-Enabled Build Completed Successfully!${NC}"
echo -e "${GREEN}======================================${NC}"
echo ""
echo -e "${BLUE}Image Information:${NC}"
docker images "${AI_IMAGE_TAG}" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
