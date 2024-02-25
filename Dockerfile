FROM python:3.12.2-slim-bullseye

WORKDIR /app
COPY requirements.txt requirements.txt

RUN pip install -r requirements.txt

COPY src run.sh ./

CMD ["sh", "run.sh"]

