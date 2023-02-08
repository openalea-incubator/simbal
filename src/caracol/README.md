# Code organisation

## Parameters

Create a file to read/write the parameters needed for the simulation.
With default values and unit (see walter?).

## Architecture

Define the global representation of several plants (a list) with an above and below ground part.
Each part is represented as an MTG. But may have different scale and properties for the above and below ground part.
As well as for between plants.

## Allocation CN

Manage the allocation (shoot & root) of C and N.

## Growth

Define the growth of the plants considering C and N allocation.

## Geometry

Reconstruct in 3D the plants' geometry 

## Environment

Compute light interception and soil N content

## Senescence

Algorithm