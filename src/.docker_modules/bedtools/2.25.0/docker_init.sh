#!/bin/sh
docker pull lbmc/bedtools:2.25.0
# docker build src/.docker_modules/bedtools/2.25.0 -t 'lbmc/bedtools:2.25.0'
# docker push lbmc/bedtools:2.25.0
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/bedtools:2.25.0" --push src/.docker_modules/bedtools/2.25.0
