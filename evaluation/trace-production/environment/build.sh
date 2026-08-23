#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
# shellcheck disable=SC1091
source "$ROOT/evaluation/trace-production/environment/versions.env"

: "${MINIMAP2_SHA256:?set MINIMAP2_SHA256 to the independently verified checksum for the official minimap2 ${MINIMAP2_VERSION} archive}"
: "${BLAST_SHA256:?set BLAST_SHA256 to the independently verified checksum for the official NCBI BLAST+ ${BLAST_VERSION} archive}"
: "${BLAST_MD5:?set BLAST_MD5 to the official NCBI checksum for the BLAST+ ${BLAST_VERSION} archive}"
: "${LEXICMAP_SHA256:?set LEXICMAP_SHA256 to the independently verified checksum for the official LexicMap ${LEXICMAP_VERSION} archive}"

docker build \
  --build-arg "MINIMAP2_VERSION=$MINIMAP2_VERSION" \
  --build-arg "MINIMAP2_TAG=$MINIMAP2_TAG" \
  --build-arg "MINIMAP2_COMMIT=$MINIMAP2_COMMIT" \
  --build-arg "MINIMAP2_SOURCE_URL=$MINIMAP2_SOURCE_URL" \
  --build-arg "MINIMAP2_SHA256=$MINIMAP2_SHA256" \
  --build-arg "BLAST_VERSION=$BLAST_VERSION" \
  --build-arg "BLAST_SOURCE_URL=$BLAST_SOURCE_URL" \
  --build-arg "BLAST_SHA256=$BLAST_SHA256" \
  --build-arg "BLAST_MD5=$BLAST_MD5" \
  --build-arg "LEXICMAP_VERSION=$LEXICMAP_VERSION" \
  --build-arg "LEXICMAP_TAG=$LEXICMAP_TAG" \
  --build-arg "LEXICMAP_COMMIT=$LEXICMAP_COMMIT" \
  --build-arg "LEXICMAP_SOURCE_URL=$LEXICMAP_SOURCE_URL" \
  --build-arg "LEXICMAP_SHA256=$LEXICMAP_SHA256" \
  -f "$ROOT/evaluation/trace-production/environment/Dockerfile" \
  -t "${COMPARATOR_IMAGE:-jam-rs-trace-comparators:$MINIMAP2_VERSION-blast-$BLAST_VERSION-lexicmap-$LEXICMAP_VERSION}" \
  "$ROOT"
