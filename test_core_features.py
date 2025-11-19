# Badge: ![Tests](https://github.com/<your_username>/pharmasight-platform/actions/workflows/test-suite.yml/badge.svg)

name: Run PharmaSight Async Test Suite

on:
  push:
  pull_request:

jobs:
  async_test_suite:
    name: Run PharmaSight Async Test Suite
    runs-on: ubuntu-latest

    steps:
      - name: Checkout repository
        uses: actions/checkout@v3

      - name: Set up Python 3.11
        uses: actions/setup-python@v4
        with:
          python-version: 3.11

      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install httpx pytest pytest-asyncio

      - name: Start Flask app in background
        run: |
          nohup python3 app.py &

      - name: Wait for server to initialize
        run: sleep 5

      - name: Run async test suite
        run: python3 test_core_features.py