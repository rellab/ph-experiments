FROM gcc:latest

WORKDIR /workspace

COPY EMpht.c .

RUN gcc EMpht.c -lm -O2 -o /usr/local/bin/empht

WORKDIR /app

CMD ["empht"]
