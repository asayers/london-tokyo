![](https://www.asayers.com/london-tokyo_small.png)

**A map of Eurasia using an oblique sinusoidal projection where the equator is the
geodesic that passes through London and Tokyo.**

[High-res version][high-res] (82MiB)

[Satellite version][satellite] (72MiB)

This map has the special property that if you draw a straight line from London
to Tokyo, you get the shortest route between those cities.
This isn't normally the case: think of those flight maps where the aeroplane
appears to follow a curved trajectory.
An aeroplane flying direct from London to Tokyo would fly straight on this map.

It uses an "oblique" projection, which means that north is not necessarily up.
In fact, the north pole is on this map, just to the right of Greenland.
You've probably already noticed that most places have the "wrong" orientation.

It's an "equal-area" map, which means everything is the correct _size_ relative
to everything else.
This is great, but the price of equal-area is that the _shapes_ of some things
must necessarily be distorted.
The question is: where to put the distortion?

The distortion is concentrated in the corners of the map:
Africa and Indonesia are squished N-S and stretched E-W.
Where it gets too bad (Tissot ratio > 2), I've masked it off in white.
That means everything you can see is within 1:2 of its actual shape.

There are two lines where things are completely undistorted:

* the horizontal line passing through London and Tokyo (via Norway and Russia)
* the vertical line running down the centre of the map (through Pakistan, Kyrgyzstan, etc.)

Click [here][tissot] (20MiB) to see some Tissot markers.

[high-res]:  https://www.asayers.com/london-tokyo.png
[satellite]: https://www.asayers.com/london-tokyo_satellite.png
[tissot]:    https://www.asayers.com/london-tokyo_tissot.png
