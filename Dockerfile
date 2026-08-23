FROM rust:1.97.1-slim-bookworm@sha256:2775a09d208ff0d7c1f50490c45b62db929e87ba1dcbc3f2132ac71a704bcdd3 AS builder

ARG JAM_RS_SOURCE_COMMIT=unknown
ENV JAM_RS_SOURCE_COMMIT=${JAM_RS_SOURCE_COMMIT}
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends make && rm -rf /var/lib/apt/lists/*
COPY Cargo.toml Cargo.lock ./
COPY src ./src
COPY benches ./benches
RUN cargo build --release --locked

FROM debian:bookworm-slim@sha256:abd67ffcfa541b485a3dff59865ab629aa048a6c613e639d36e7456b0b229241

LABEL org.opencontainers.image.source="https://github.com/St4NNi/jam-rs"
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl && rm -rf /var/lib/apt/lists/*
COPY --from=builder /app/target/release/jam /usr/local/bin/jam
WORKDIR /work
USER 65532:65532
ENTRYPOINT ["/usr/local/bin/jam"]
