P0-T4 answer to open question **Q1**: what state is ``GSE200563_processed_data.txt.gz``
in? Left panel is the decisive test — Q3 normalisation is precisely the operation
that equalises each AOI's third quartile, so a near-zero coefficient of variation
in column-wise Q3 is its signature. The middle panel shows library size over the
same AOIs for contrast: it varies normally, which is what rules out the
possibility that the flat Q3 is an artefact of uniform input. The right panel
places the value distribution, separating a raw-count scale from a log scale.

``NegProbe-WTX`` is excluded from these statistics; it is a control, not a gene.

The verdict drives ``normalisation.method`` in ``config/config.yaml``: if the
matrix arrives already Q3-normalised, P0-T5 *verifies* rather than *applies*.
