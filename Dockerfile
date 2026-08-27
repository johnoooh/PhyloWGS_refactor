FROM --platform=$BUILDPLATFORM golang:1.22-bookworm AS build

ARG TARGETOS
ARG TARGETARCH

WORKDIR /src
COPY . .
RUN CGO_ENABLED=1 GOOS=${TARGETOS} GOARCH=${TARGETARCH} go build -o /out/phylowgs-go .

FROM debian:bookworm-slim

ARG PHYLOWGS_GO_TAG=go-port
ENV PHYLOWGS_GO_TAG=${PHYLOWGS_GO_TAG}

LABEL org.opencontainers.image.title="phylowgs-go" \
      org.opencontainers.image.description="Pure Go reimplementation of the PhyloWGS MCMC sampler (go-port branch)." \
      org.opencontainers.image.source="https://github.com/mskcc/phylowgs" \
      org.opencontainers.image.authors="John Orgera (orgeraj@mskcc.com)" \
      org.opencontainers.image.version=${PHYLOWGS_GO_TAG}

RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY --from=build /out/phylowgs-go /usr/local/bin/phylowgs-go
