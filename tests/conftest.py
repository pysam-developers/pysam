import os


def pytest_report_header(config):
    text = []

    if "REF_PATH" in os.environ:
        text.append("pysam: overriding REF_PATH to disable external reference lookups")
    os.environ["REF_PATH"] = ":"

    return text
