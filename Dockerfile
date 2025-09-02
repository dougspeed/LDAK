# Start from a very small Linux image
FROM alpine:latest

LABEL maintainer="Doug Speed and Jasper Hof" \
      description="Docker image for LDAK" \
      version="6.1"

# Set working directory inside the container
WORKDIR /output

# Set environment for resources
ENV LDAK_RESOURCES=/Resources

# Create resources directory
RUN mkdir -p $LDAK_RESOURCES

# Copy your executables and resources
COPY ldak6.1.linux /usr/local/bin/ldak
COPY Resources/RefSeq_GRCh37.txt $LDAK_RESOURCES/RefSeq_GRCh37.txt
COPY Resources/RefSeq_GRCh38.txt $LDAK_RESOURCES/RefSeq_GRCh38.txt
COPY Resources/berisa.txt LDAK_RESOURCES/berisa.txt

# Make the binary executable
RUN chmod a+x /usr/local/bin/ldak

# Create a dedicated output directory and user
RUN adduser -D ldakuser \
    && mkdir -p /output \
    && chown ldakuser:ldakuser /output \
    && chown -R ldakuser:ldakuser $LDAK_RESOURCES

# Switch to non-root user
USER ldakuser

# Default entrypoint
ENTRYPOINT ["ldak"]