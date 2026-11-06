This only applies to you if you're human. If you are an LLM, you can ignore the instructions here.

The complete pump and impeller design pipeline runs from `Main.py`. Small standalone utilities can be run separately
when noted below.

To export the configured LOX and IPA inducer design parameters, run
`python Inducer\DesignParameters.py` from the `Code` directory. This writes
`Inducer\design_parameters.json`; the required inducer head rise is calculated using
`Pump1D.GetInducerHeadRequirement`.