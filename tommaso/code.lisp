
(defun hello-world ()
  (print "Hello world"))


(defun vstrands (j)
  (cdr (assoc :vstrands j)))

(defun strands (vstrands i)
  (elt vstrands i))


(defun length-of-strand (strands)
  (length (cdr (assoc :stap strands))))




#|
(length-of-strand (strands (vstrands *j*) 1))

(apply #'+ (mapcar #'length-of-strand (vstrands *j*)))

(vstrands *j*)

(assoc :stap (vstrands *j*))

|#
