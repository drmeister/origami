
(defun LIST-OF-VSTRANDS-FROM-JSON (json)
  (cdr (assoc "vstrands" json)))


(defclass vstrand ()
  ((staple-json :initarg :staple-json :reader staple-json)
   (scaffold-json :initarg :scaffold-json :reader scaffold-json)
   ...
   (staple-vec :initarg :staple-vec :reader staple-vec)
   (scaffold-vec :initarg :scaffold-vec :reader scaffold-vec)
   ...))

   
(defun parse-json (json)
  (let ((vstrands (make-hash-table :test #'eql)))
    (loop for json-vstrand-info in (LIST-OF-VSTRANDS-FROM-JSON json)
       for num = (cdr (assoc "num" json-vstrand-info))
       for col = (cdr (assoc "col" json-vstrand-info))
       for row = (cdr (assoc "row" json-vstrand-info))
       for staple-json = (cdr (assoc "stap" json-vstrand-info))
       for scaffold-json = (cdr (assoc "scaf" json-vstrand-info))
       for skip = (cdr (assoc "skip" json-vstrand-info))
       for loop = (cdr (assoc "loop" json-vstrand-info))
       for scaf-loop = (cdr (assoc "scafLoop" json-vstrand-info))
       for stap-loop = (cdr (assoc "stapLoop" json-vstrand-info))
       for vstrand = (progn
		       (or (= (length staple) (length scaffold))
			   (error "In vstrand ~a the staple length ~a and the scaffold length ~a are not the same - and they must be"
				  num
				  (length staple)
				  (length scaffold)))
		       (make-instance 'vstrand
				      :staple-json staple-json
				      :scaffold-json scaffold-json
				      :skip-json skip
				      :loop-json loop
				      :scaf-loop-json scaf-loop
				      :stap-loop-json stap-loop
				      :staple-vec (make-array (length staple-json))
				      :scaffold-vec (make-array (length scaffold-json))
				      :row row
				      :column col))
       do (setf (gethash num vstrands) vstrand)
       do (BUILD-PRE-NODE (staple-vec vstrand) staple-json)
       do (BUILD-PRE-NODE (scaffold-vec vstrand) scaffold-json)
       do (INTRA-HELIX-CONNECT (staple-vec vstrand) (scaffold-vec vstrand)))
    ;; Connect the prenodes
    (loop for vstrand being the hash-values in vstrands using (hash-key num)
       for staple-json = (staple-json vstrand)
       for staple-vec = (staple-vec vstrand)
       for scaffold-json = (scaffold-json vstrand)
       for scaffold-vec = (scaffold-vec vstrand)
       do (CONNECT-EVERYTHING vstrands num staple-json staple-vec #'staple-vec)
       do (CONNECT-EVERYTHING vstrands num scaffold-json scaffold-vec #'scaffold-vec))
    ;; Deal with loop



    
       (defun CONNECT-EVERYTHING (vstrands num json vec accessor)
	 (loop for json-info in json
	    for index from 0 below (length json)
	    for node = (elt vec index)
	    when node
	    for forward-node = (LOOKUP-NODE vstrands (first json-info) (second json-info) accessor)
	    for back-node = (LOOKUP-NODE vstrands (third json-info) (fourth json-info) accessor)
	    do (when forward-node (setf (forward-node node) forward-node))
	    do (when back-node (setf (back-node node) back-node))))
	      
       (defclass pre-node ()
	 ((hbond-node :initarg :hbond-node :accessor hbond-node)
	  (forward-node :initarg :forward-node :accessor forward-node)
	  (back-node :initarg :back-node :accessor back-node)))
       


       (defun LOOKUP-NODE (vstrands vstrand-num pos accessor)
	 (let* ((vstrand (gethash vstrand-num vstrands))
		(strand (funcall accessor vstrand)))
	   (elt strand pos)))

       
    vstrands))


	 
				    
	 collect row)
	 
    
