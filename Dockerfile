FROM rust:1-slim AS builder
WORKDIR /app
COPY . .
RUN cargo build --release --locked

FROM debian:stable-slim
COPY --from=builder /app/target/release/jam /usr/local/bin/jam
ENTRYPOINT ["jam"]
