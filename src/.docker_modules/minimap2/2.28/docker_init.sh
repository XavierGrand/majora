#!/bin/sh
# docker pull xgrand/minimap2:2.28
docker build src/.docker_modules/minimap2/2.28 -t 'xgrand/minimap2:2.28'
docker push xgrand/minimap2:2.28
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/minimap2:2.28" --push src/.docker_modules/minimap2/2.28
