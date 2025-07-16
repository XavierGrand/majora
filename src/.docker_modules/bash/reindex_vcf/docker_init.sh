#!/bin/sh
docker pull xgrand/reindex_vcf:1.1
# docker build src/.docker_modules/bash/reindex_vcf -t 'xgrand/reindex_vcf:1.1'
# docker push xgrand/reindex_vcf:1.1
docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/reindex_vcf:1.1" --push src/.docker_modules/bash/reindex_vcf 
