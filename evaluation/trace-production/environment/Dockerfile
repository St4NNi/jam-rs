# syntax=docker/dockerfile:1.7
FROM debian:bookworm-slim AS build

ARG MINIMAP2_VERSION=2.31
ARG MINIMAP2_TAG=v2.31
ARG MINIMAP2_COMMIT=3c28777
ARG MINIMAP2_SOURCE_URL=https://github.com/lh3/minimap2/releases/download/v2.31/minimap2-2.31_x64-linux.tar.bz2
ARG MINIMAP2_SHA256
ARG BLAST_VERSION=2.17.0
ARG BLAST_SOURCE_URL=https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-x64-linux.tar.gz
ARG BLAST_SHA256
ARG BLAST_MD5
ARG LEXICMAP_VERSION=0.9.0
ARG LEXICMAP_TAG=v0.9.0
ARG LEXICMAP_COMMIT=e8dec8d
ARG LEXICMAP_SOURCE_URL=https://github.com/shenwei356/LexicMap/releases/download/v0.9.0/lexicmap_linux_amd64.tar.gz
ARG LEXICMAP_SHA256

RUN test -n "$MINIMAP2_SHA256" || (echo "MINIMAP2_SHA256 is required and must be independently verified" >&2; exit 64) \
    && test -n "$BLAST_SHA256" || (echo "BLAST_SHA256 is required and must be independently verified" >&2; exit 64) \
    && test -n "$BLAST_MD5" || (echo "BLAST_MD5 is required and must be independently verified" >&2; exit 64) \
    && test -n "$LEXICMAP_SHA256" || (echo "LEXICMAP_SHA256 is required and must be independently verified" >&2; exit 64)

RUN apt-get update \
    && apt-get install --no-install-recommends -y ca-certificates curl bzip2 zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
RUN curl --fail --location --silent --show-error "$MINIMAP2_SOURCE_URL" --output minimap2.tar.bz2 \
    && printf '%s  %s\n' "$MINIMAP2_SHA256" minimap2.tar.bz2 | sha256sum --check --strict - \
    && mkdir minimap2 \
    && tar --extract --bzip2 --file minimap2.tar.bz2 --directory minimap2 \
    && minimap2_bin="$(find minimap2 -type f -name minimap2 -perm -u+x -print -quit)" \
    && test -n "$minimap2_bin" \
    && install -D -m 0755 "$minimap2_bin" "/opt/comparators/minimap2-$MINIMAP2_VERSION"

RUN curl --fail --location --silent --show-error "$BLAST_SOURCE_URL" --output blast.tar.gz \
    && printf '%s  %s\n' "$BLAST_SHA256" blast.tar.gz | sha256sum --check --strict - \
    && printf '%s  %s\n' "$BLAST_MD5" blast.tar.gz | md5sum --check --strict - \
    && mkdir blast \
    && tar --extract --gzip --file blast.tar.gz --directory blast --strip-components=1 \
    && install -D -m 0755 blast/bin/blastn "/opt/comparators/blastn-$BLAST_VERSION"

RUN curl --fail --location --silent --show-error "$LEXICMAP_SOURCE_URL" --output lexicmap.tar.gz \
    && printf '%s  %s\n' "$LEXICMAP_SHA256" lexicmap.tar.gz | sha256sum --check --strict - \
    && mkdir lexicmap \
    && tar --extract --gzip --file lexicmap.tar.gz --directory lexicmap \
    && lexicmap_bin="$(find lexicmap -type f -name lexicmap -perm -u+x -print -quit)" \
    && test -n "$lexicmap_bin" \
    && install -D -m 0755 "$lexicmap_bin" "/opt/comparators/lexicmap-$LEXICMAP_VERSION"

FROM debian:bookworm-slim

ARG MINIMAP2_VERSION=2.31
ARG MINIMAP2_TAG=v2.31
ARG MINIMAP2_COMMIT=3c28777
ARG BLAST_VERSION=2.17.0
ARG LEXICMAP_VERSION=0.9.0
ARG LEXICMAP_TAG=v0.9.0
ARG LEXICMAP_COMMIT=e8dec8d

LABEL org.opencontainers.image.description="Pinned BLAST+, minimap2, and LexicMap comparator environment" \
      org.jam-rs.minimap2.version="$MINIMAP2_VERSION" \
      org.jam-rs.minimap2.tag="$MINIMAP2_TAG" \
      org.jam-rs.minimap2.commit="$MINIMAP2_COMMIT" \
      org.jam-rs.blast.version="$BLAST_VERSION" \
      org.jam-rs.lexicmap.version="$LEXICMAP_VERSION" \
      org.jam-rs.lexicmap.tag="$LEXICMAP_TAG" \
      org.jam-rs.lexicmap.commit="$LEXICMAP_COMMIT"

RUN apt-get update \
    && apt-get install --no-install-recommends -y ca-certificates libgomp1 zlib1g \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd --system --gid 10001 comparator \
    && useradd --system --uid 10001 --gid comparator --create-home --home-dir /home/comparator comparator

COPY --from=build "/opt/comparators/minimap2-$MINIMAP2_VERSION" /usr/local/bin/minimap2
COPY --from=build "/opt/comparators/blastn-$BLAST_VERSION" /usr/local/bin/blastn
COPY --from=build "/opt/comparators/lexicmap-$LEXICMAP_VERSION" /usr/local/bin/lexicmap

USER comparator
WORKDIR /work
ENTRYPOINT ["/bin/sh"]
