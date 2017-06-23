
(in-package :cl-jupyter-user)

(defparameter *s* #P"cadnano_test_file/super_barcode_hex.json")


(defparameter *s* #P"cadnano_test_file/gap_vs_skip.json")
(defparameter *origami* (parse-cadnano *s*))

*origami*
(fill-nodes-with-residues *origami*)

*origami*

(defparameter *m* (build-one-molecule *origami*))

(chem:save-mol2 *m* "test.mol2")
