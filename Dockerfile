FROM python:3.12.2-slim-bullseye

WORKDIR /app
COPY requirements.txt requirements.txt

RUN pip install -r requirements.txt

COPY src ./

ENV PORT=8080

CMD uvicorn app:app --host 0.0.0.0 --port $PORT
