# ExaHyPE 2 Training

## Installation

Please install [Docker](https://docs.docker.com/engine/install/), launch Docker and then run

```bash
docker pull peanoframework/training
```

## Training

After installation, run

```bash
docker run --rm -p 9999:9999 peanoframework/training
```

You should see

```bash
http://127.0.0.1:9999/lab?token=
```

Click on that link (ctrl+left-click) or enter the link in the address bar of your web browser.

Then use the navigation bar to open the exercises (e.g., [acoustic](acoustic)).

## Interactive Docker

You can also run an interactive Docker container and mount the training repository:

```bash
docker run -it -v ${PWD}:/training --rm -p 9999:9999 peanoframework/training
```

Then you may start the jupyter lab session in the `/training` folder:

```bash
jupyter lab --allow-root --port=9999 --no-browser --ip=0.0.0.0
```
