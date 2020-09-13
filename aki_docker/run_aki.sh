#!/bin/bash
docker run --rm --name 4ce -d -v /mnt/c/Users/Mummy/Documents/4CE/aki_docker/data:/4ceData -p 8787:8787 -p 2200:22 -e CONTAINER_USER_USERNAME=byorntan -e CONTAINER_USER_PASSWORD=rushhour2020 dbmi/4ce-analysis:v1.0aki

