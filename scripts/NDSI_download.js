//VERSION=3
//Subset of S2L2A raw bands for calculating NDSI, original data (no harmonization)

function setup() {
  return {
    input: [{
      bands: ["B03", "B11", "B04", "B02"],
      units: "DN"
    }],
    output: {
      id: "default",
      bands: 4,
      sampleType: SampleType.UINT16
    }
  }
}

function evaluatePixel(sample) {
    return [ sample.B02, sample.B03, sample.B04, sample.B11]
}
