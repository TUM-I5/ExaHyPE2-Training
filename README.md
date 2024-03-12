# ExaHyPE 2 Training

## Download the Instructions

The PDF is uploaded from the CI pipeline as an artifact.

Go to the `Actions` tab, click on the latest pipeline run and download the `ExaHyPE2_Training` artifact.

Extract it and open the PDF.

## Installation

Please install [Docker](https://docs.docker.com/engine/install/), launch Docker and then run

```bash
docker pull peanoframework/training
```

## Training

After installation, clone this repository

```bash
git clone https://github.com/TUM-I5/ExaHyPE2-Training.git
```

Navigate into the cloned repository and run

```bash
docker run -v ${PWD}:/training --rm -p 9999:9999 peanoframework/training
```

You should see

```bash
http://127.0.0.1:9999/lab?token=
```

Click on that link (ctrl+left-click) or enter the link in the address bar of your web browser.

Then use the navigation bar to open the exercises (e.g., [acoustic](acoustic)).
