
(setf *default-pathname-defaults* #P"~/cando/origami/tommaso/")

(load-*bases* #P"~/cando/origami/tommaso/base-pair-pdb-file/")

;;;(defparameter *s* #P"~/cando/origami/tommaso/cadnano_test_file/super_barcode_hex.json")

(defparameter *s* #P"~/cando/origami/tommaso/cadnano_test_file/gap_vs_skip.json")
(defparameter *origami* (parse-cadnano *s*))

*origami*
(fill-nodes-with-residues *origami*)

*origami*

(defparameter *m* (build-one-molecule *origami*))

(chem:save-mol2 *m* (namestring (merge-pathnames #P"test.mol2")))


(defparameter *staple* (build-energy-rigid-body-staple *origami*))

