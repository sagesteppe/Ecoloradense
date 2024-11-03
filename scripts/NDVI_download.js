//VERSION=3
//Subset of S2L2A raw bands for calculating NDSI, original data (no harmonization)

function setup() {
  return {
    input: [{
      bands: ["B04", "B08"],
      units: "DN"
    }],
    output: {
      id: "default",
      bands: 2,
      sampleType: SampleType.UINT16
    }
  }
}

function evaluatePixel(sample) {
    return [sample.B04, sample.B08]
}
