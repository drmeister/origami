
(defun vstrands (j)
  (cdr (assoc :vstrands j)))

(defun strands (vstrands i)
  (elt vstrands i))

(defun length-of-strand (strands)
  (length (cdr (assoc :stap strands))))

#|
(length-of-strand (strands (vstrands *j*) 1))

(apply #'+ (mapcar #'length-of-strand (vstrands *j*)))

(assoc :stap (vstrands *j*))

|#

(defclass pre-node () (ID ;must the node be identified with a numerical label?
			fwd-strand
			fwr-pos
			bkd-strand
			bkd-pos
			fwd-node
			bkw-node
			hb-node))
(defclass pre-spec-node () (ID :ID)
			(fwd-strand :fwd-strand)
			(fwr-pos :fwr-pos)
			(bkd-strand :bkd-strand)
			(bkd-pos :bkd-pos)
			(fwd-node :fwd-node)
			(bkw-node :bkw-node)
			(hb-node :hb-node)))


(defclass node () (ID
			fwd-node
			bkw-node
			hb-node))

(defun make-strand-of-node (lengthST)
  (make-array lenghtST :initial-element nil :element-type 'pre-node))

;(make-instance pre-node)
(defun fill-node (pnode coordinate content)
  (setf (slot-value pnode 'coordinate) content))

;(defun fill-spec-node (pnode coordinate content)
(make-instance  content)
;)


