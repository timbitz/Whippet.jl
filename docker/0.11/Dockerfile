##### BASE IMAGE #####
FROM ubuntu:18.04

# install dependencies
RUN apt-get update \
  && apt-get install -y tzdata \
  && ln -fs /usr/share/zoneinfo/Europe/Berlin /etc/localtime \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get install --yes make g++ libgslcblas0 libglib2.0-dev libgtk2.0 wget unzip

# install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz \
  && tar xzvf julia-0.6.4-linux-x86_64.tar.gz \
  && ln -s /julia-9d11f62bcb/bin/julia /usr/bin/julia 

# install Whippet.jl
RUN julia -e 'Pkg.add("Whippet")' \
  && julia -e 'using Whippet'

# Clean up
RUN rm /julia-0.6.4-linux-x86_64.tar.gz

WORKDIR /root/.julia/v0.6/Whippet/bin/
