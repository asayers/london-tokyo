#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
from PIL import Image


# ---------------------------------------------------------------------------
# Configuration

LONDON = (51.5074, -0.1278)       # latitude, longitude
TOKYO  = (35.6762, 139.6503)

# Approximate centre of the Bering Strait.  The top edge of the output passes
# through this point.
BERING_STRAIT = (65.75, -168.97694)

# Horizontal projected margin beyond London and Tokyo.
MARGIN = 22.0

# Bottom edge in oblique latitude.
YMIN = -65.0

# Tissot major/minor axis ratio:
#   <= FADE_START   unchanged
#   >= FADE_END     white
#   between         smooth fade to white
FADE_START = 1.5
FADE_END = 2.5
# FADE_START = 1000.0
# FADE_END = 1000.0
# FADE_START = 2.0
# FADE_END = 3.0

# Width of rendered image. Height follows from the projection bounds.
OUTPUT_WIDTH = 4096


# ---------------------------------------------------------------------------
# Spherical geometry

def unit_vector(lat, lon):
    """Geographic latitude/longitude -> Earth-centred unit vector."""

    lat, lon = np.radians([lat, lon])

    return np.array([
        np.cos(lat) * np.cos(lon),
        np.cos(lat) * np.sin(lon),
        np.sin(lat),
    ])


# London and Tokyo as Earth-centred unit vectors.
L = unit_vector(*LONDON)
T = unit_vector(*TOKYO)

# N is the pole of the London–Tokyo great circle.
#
# Therefore that great circle becomes phi'=0 in our oblique coordinate
# system.
N = np.cross(L, T)
N /= np.linalg.norm(N)

# Initially choose London as lambda'=0.
X = L

# Complete the right-handed orthonormal basis.
Y = np.cross(N, X)
Y /= np.linalg.norm(Y)


def geographic_to_oblique(lat, lon):
    """
    Geographic lat/lon -> oblique spherical lambda/phi, in degrees.

    phi=0 is the London–Tokyo great circle.
    lambda=0 is London.
    """

    p = unit_vector(lat, lon)

    lam = np.arctan2(
        p @ Y,
        p @ X,
    )

    phi = np.arcsin(
        np.clip(p @ N, -1.0, 1.0)
    )

    return np.degrees(lam), np.degrees(phi)


# Tokyo's lambda is the angular length of the London–Tokyo great-circle arc.
TOKYO_LAMBDA, _ = geographic_to_oblique(*TOKYO)

# Put lambda=0 for the actual map halfway between London and Tokyo.
#
# London and Tokyo are consequently at equal and opposite x coordinates.
LAMBDA_MID = TOKYO_LAMBDA / 2.0

LONDON_X = -LAMBDA_MID
TOKYO_X = +LAMBDA_MID


# ---------------------------------------------------------------------------
# Projection
#
# Midpoint-centred oblique sinusoidal:
#
#     mu = lambda' - lambda_mid
#
#     x = mu cos(phi')
#     y = phi'
#
# This projection is equal-area everywhere.
#
# The London–Tokyo geodesic is y=0.  Along it the projection has unit scale
# and zero angular distortion.
#
# The perpendicular great circle through the London–Tokyo midpoint is x=0
# and likewise has unit scale.


def inverse_projection(x, y):
    """
    Projected x/y -> geographic latitude/longitude.

    x and y use degree-equivalent units.
    """

    phi = np.radians(y)

    cos_phi = np.cos(phi)

    # Invert:
    #
    #     x = mu cos(phi)
    #
    mu = x / cos_phi

    # Restore the original oblique longitude, whose zero was London.
    lam = np.radians(mu + LAMBDA_MID)

    # Oblique spherical coordinates -> Earth-centred vector.
    p = (
        np.cos(phi)[..., None]
        * np.cos(lam)[..., None]
        * X

        +

        np.cos(phi)[..., None]
        * np.sin(lam)[..., None]
        * Y

        +

        np.sin(phi)[..., None]
        * N
    )

    lat = np.degrees(
        np.arcsin(
            np.clip(p[..., 2], -1.0, 1.0)
        )
    )

    lon = np.degrees(
        np.arctan2(
            p[..., 1],
            p[..., 0],
        )
    )

    return lat, lon


# ---------------------------------------------------------------------------
# Distortion
#
# In an orthonormal tangent basis on the sphere, the differential of the
# sinusoidal projection is
#
#           [ 1   -mu sin(phi) ]
#     J  =  [                  ]
#           [ 0        1       ]
#
# where mu is measured in radians.
#
# det(J)=1, expressing the equal-area property.
#
# The singular values of J are the major/minor radii of the local Tissot
# ellipse. Since their product is one, their ratio completely describes
# the local shape distortion.


def tissot_axis_ratio(x, y):
    """
    Return major/minor axis ratio of the local Tissot indicatrix.

    1.0 = locally conformal
    2.0 = infinitesimal circle becomes a 2:1 ellipse
    etc.
    """

    phi = np.radians(y)

    # mu must be radians here.
    mu = np.radians(
        x / np.cos(phi)
    )

    q = np.abs(
        mu * np.sin(phi)
    )

    # Largest singular value of
    #
    #     [[1, q],
    #      [0, 1]]
    #
    smax = (
        np.sqrt(q * q + 4.0) + q
    ) / 2.0

    # det(J)=1, so:
    #
    #     smin = 1/smax
    #
    # hence the major/minor ratio is smax².
    return smax * smax


# ---------------------------------------------------------------------------
# Bilinear raster sampling

def sample_equirectangular(source, lat, lon):
    """
    Bilinearly sample a global equirectangular raster.

    Assumes:
        horizontal: longitude -180 .. +180
        vertical:   latitude   +90 .. -90
    """

    h, w, _ = source.shape

    u = (lon + 180.0) / 360.0 * (w - 1)
    v = (90.0 - lat) / 180.0 * (h - 1)

    # Longitude wraps around.
    u %= w

    # Latitude does not.
    v = np.clip(v, 0, h - 1)

    u0 = np.floor(u).astype(np.int32)
    v0 = np.floor(v).astype(np.int32)

    u1 = (u0 + 1) % w
    v1 = np.minimum(v0 + 1, h - 1)

    du = (u - u0)[..., None]
    dv = (v - v0)[..., None]

    return (
        (1 - du) * (1 - dv) * source[v0, u0]
        + du * (1 - dv) * source[v0, u1]
        + (1 - du) * dv * source[v1, u0]
        + du * dv * source[v1, u1]
    )


# ---------------------------------------------------------------------------
# Render

def render(source_path, output_path, width=OUTPUT_WIDTH):

    # Pillow protects against accidentally enormous/decompression-bomb images.
    # Blue Marble is intentionally large.
    Image.MAX_IMAGE_PIXELS = None

    source = np.asarray(
        Image.open(source_path).convert("RGB"),
        dtype=np.float32,
    )

    # -----------------------------------------------------------------------
    # Bounds

    xmin = LONDON_X - MARGIN
    xmax = TOKYO_X + MARGIN

    # Top edge passes through the centre of the Bering Strait.
    _, ymax = geographic_to_oblique(*BERING_STRAIT)

    ymin = YMIN

    # x and y use the same degree-equivalent map units, so preserve that
    # aspect ratio in the raster dimensions.
    height = round(
        width
        * (ymax - ymin)
        / (xmax - xmin)
    )

    print(f"London/Tokyo separation: {TOKYO_LAMBDA:.3f}°")
    print(f"London x: {LONDON_X:.3f}°")
    print(f"Tokyo x:  {TOKYO_X:.3f}°")
    print(f"Top edge: {ymax:.3f}°")
    print(f"Bottom:   {ymin:.3f}°")
    print(f"Output:   {width} × {height}")

    # -----------------------------------------------------------------------
    # Output grid

    xs = np.linspace(
        xmin,
        xmax,
        width,
        dtype=np.float64,
    )

    # Raster rows run top -> bottom.
    ys = np.linspace(
        ymax,
        ymin,
        height,
        dtype=np.float64,
    )

    xx, yy = np.meshgrid(xs, ys)

    # -----------------------------------------------------------------------
    # Inverse projection

    lat, lon = inverse_projection(xx, yy)

    # -----------------------------------------------------------------------
    # Sample source raster

    image = sample_equirectangular(
        source,
        lat,
        lon,
    )

    # -----------------------------------------------------------------------
    # Distortion fade

    distortion = tissot_axis_ratio(xx, yy)

    # Normalised fade:
    #
    #   distortion <= 2  -> 0
    #   distortion >= 3  -> 1
    #
    t = np.clip(
        (distortion - FADE_START)
        / (FADE_END - FADE_START),
        0.0,
        1.0,
    )

    # Smoothstep:
    #
    # smoother than a linear alpha ramp, with zero derivative at each end.
    t = t * t * (3.0 - 2.0 * t)

    # Composite directly onto white rather than using transparency.
    #
    # This means the PNG displays identically regardless of viewer
    # background.
    image = (
        image * (1.0 - t[..., None])
        + 255.0 * t[..., None]
    )

    image = np.clip(
        image,
        0,
        255,
    ).astype(np.uint8)

    # -----------------------------------------------------------------------
    # Save directly with Pillow.
    #
    # No matplotlib means:
    #   - no border
    #   - no title
    #   - no axes
    #   - no padding
    #   - exact raster dimensions

    Image.fromarray(
        image,
        mode="RGB",
    ).save(
        output_path,
        quality=95,
    )


# ---------------------------------------------------------------------------
# CLI

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Render an equal-area map centred on the "
            "London–Tokyo great-circle route."
        )
    )

    parser.add_argument(
        "source",
        type=Path,
        help=(
            "Global equirectangular source raster "
            "(e.g. NASA Blue Marble)."
        ),
    )

    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        default=Path("london-tokyo.png"),
    )

    parser.add_argument(
        "--width",
        type=int,
        default=OUTPUT_WIDTH,
        help=f"Output width in pixels (default {OUTPUT_WIDTH}).",
    )

    args = parser.parse_args()

    render(
        args.source,
        args.output,
        args.width,
    )


if __name__ == "__main__":
    main()
