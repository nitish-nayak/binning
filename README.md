playing around with some binning optimization algorithms (requires c++17)

- `interfaces/root.h` has some applications with ROOT (requires ROOT compiled with c++17, ideally `6.30+`)

### Minimal setup for root
```bash
# use a container with a modern root version and mount our repo to it
docker pull rootproject/root
docker run --rm -it -v $(pwd):/binning rootproject/root bash

# then inside container
cd binning/interfaces/tests
make roottest
./roottest
```
