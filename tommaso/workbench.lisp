
(ql:quickload "cl-json")
(ql:quickload "read-csv")

(setf *default-pathname-defaults* #P"~/cando/origami/tommaso/")

(load (compile-file "parser-cadnano.lisp"))

(load-*bases* #P"~/cando/origami/tommaso/base-pair-pdb-file/")

(defparameter *s* #P"~/cando/origami/tommaso/cadnano_test_file/super_barcode_hex.json")

;;;(defparameter *s* #P"~/cando/origami/tommaso/cadnano_test_file/dna-2-nicks.json")
(defparameter *origami* (parse-cadnano *s*))
(fill-nodes-with-residues *origami*)



(defparameter *staple* (build-energy-rigid-body-staple *origami*))
(defparameter *nonbond* (build-energy-rigid-body-nonbond *origami*))

(defparameter *ef* (build-rigid-body-energy-function *origami*))
(chem:rigid-body-energy-function-add-term *ef* *staple*)


(chem:rigid-body-energy-function-add-term *ef* *nonbond*)

(chem:dump-terms *ef*)
(format t "Done~%")
(defparameter *min* (chem:make-minimizer *ef*))

(format t "~%quat-trans at 0 -> ~a~%" (chem:rigid-body-energy-function-get-position *ef* 1))
(format t "Done setup~%")

(chem:enable-print-intermediate-results *min*)

(randomize-rigid-body-energy-function-position *ef*)
(let ((mol (build-one-aggregate *origami* *ef*)))
  (cando:save-mol2 mol "before.mol2"))

(chem:enable-debug *staple*)
(chem:enable-debug *nonbond*)
(chem:disable-debug *staple*)
(chem:disable-debug *nonbond*)

(chem:calculate-energy-and-force *ef*)

(progn
  (chem:set-maximum-number-of-steepest-descent-steps *min* 4000)
  (chem:set-maximum-number-of-truncated-newton-steps *min* 0)
  (chem:set-steepest-descent-tolerance *min* 500.0)
  (chem:set-conjugate-gradient-tolerance *min* 0.01))

(chem:minimize *min*)


(let ((mol (build-one-aggregate *origami* *ef*)))
  (cando:save-mol2 mol "after.mol2"))
