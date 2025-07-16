#!/bin/sh
# docker pull xgrand/dorado:0.8.2
docker build src/.docker_modules/dorado/0.8.2 -t 'xgrand/dorado:0.8.2'
docker push xgrand/dorado:0.8.2
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/dorado:0.8.2" --push src/.docker_modules/dorado/0.8.2