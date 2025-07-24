FROM biopipen/base:dev

# in case we have update dependencies
COPY tests/conda/env_biopipen.yml /biopipen/tests/conda/env_biopipen.yml

# Install dependencies
RUN micromamba env update -n base -f /biopipen/tests/conda/env_biopipen.yml && \
    micromamba clean --all --yes

# Test it by:
# docker run -it --rm biopipen/test:dev \
#   pipen run vcf Vcf2Bed --in.infile /biopipen/tests/test_vcf/Vcf2Bed/data/tobed.vcf
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh" ]
