(in-package :origami)

(defun LIST-OF-VSTRANDS-FROM-JSON (json)
  (cdr (assoc :vstrands json)))


(defclass vstrand ()
  ((staple-json :initarg :staple-json :reader staple-json)
   (scaffold-json :initarg :scaffold-json :reader scaffold-json)
   (skip-json :initarg :skip-json :reader skip-json)
   (loop-json :initarg :loop-json :reader loop-json)
   (scaf-loop-json :initarg :scaf-loop-json :reader scaf-loop-json)
   (stap-loop-json :initarg :stap-loop-json :reader stap-loop-json)
   (staple-vec :initarg :staple-vec :reader staple-vec)
   (scaffold-vec :initarg :scaffold-vec :reader scaffold-vec)
   (p-strand :initarg :p-strand :accessor p-strand)
   (n-strand :initarg :n-strand :accessor n-strand)
   (row :initarg :row :reader row)
   (col :initarg :column :reader column)
   ))

  ;check if not ("format":"3.0")
(defun parse-json (json &optional csv-data)
  (format t "Starting parse-json~%")
  (let ((vstrands (make-hash-table :test #'eql)))
    (loop for json-vstrand-info in (LIST-OF-VSTRANDS-FROM-JSON json)
       for num = (cdr (assoc :num json-vstrand-info))
       for col = (cdr (assoc :col json-vstrand-info))
       for row = (cdr (assoc :row json-vstrand-info))
       for staple-json = (cdr (assoc :stap json-vstrand-info))
       for scaffold-json = (cdr (assoc :scaf json-vstrand-info))
       for skip = (cdr (assoc :skip json-vstrand-info))
       for loop = (cdr (assoc :loop json-vstrand-info))
       for scaf-loop = (cdr (assoc :scafLoop json-vstrand-info))
       for stap-loop = (cdr (assoc :stapLoop json-vstrand-info))
       for vstrand = (progn
                       (or (= (length staple-json) (length scaffold-json))
                           (error "In vstrand ~a the staple length ~a and the scaffold length ~a are not the same - and they must be"
                                  num
                                  (length staple-json)
                                  (length scaffold-json)))
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
       do (format t "vstrand -> ~a~%" vstrand)
       do (setf (gethash num vstrands) vstrand)
       do (BUILD-NODE staple-json (staple-vec vstrand) num :sT)
       do (BUILD-NODE scaffold-json (scaffold-vec vstrand) num :sC)
       do (intra-helix-connect (staple-vec vstrand) (scaffold-vec vstrand))
         )
    ;; connect the nodes
    (loop for vstrand being the hash-values in vstrands using (hash-key num)
       for staple-json = (staple-json vstrand)
       for staple-vec = (staple-vec vstrand)
       for scaffold-json = (scaffold-json vstrand)
       for scaffold-vec = (scaffold-vec vstrand)
       do (CONNECT-EVERYTHING vstrands num staple-json staple-vec #'staple-vec)
       do (CONNECT-EVERYTHING vstrands num scaffold-json scaffold-vec #'scaffold-vec))
    ;;deal with sequence
    (format t "here~%")
    (unless csv-data
      (warn "No CSV; nucleobase information will be defaulted"))
    (format t "About to assign node-staple-sequences~%")
    (let ((node-staple-sequences (make-hash-table :test #'eq)))
      (loop for (num pos sequence) in csv-data
         for start-node = (lookup-node vstrands num pos 'staple-vec)
         do (setf (gethash start-node node-staple-sequences) sequence))
      ;; Deal with loop
      ;; (skip-loop skip loop staple-vec scaffold-vec) -> return new staple and scaffold
      (loop for vstrand being the hash-values in vstrands using (hash-key num)
         for loop-json = (loop-json vstrand)
         for staple-vec = (staple-vec vstrand)
         for skip-json = (skip-json vstrand)
         for scaffold-vec = (scaffold-vec vstrand)
         do (format t "duplex num ~a~%" num)
         do (multiple-value-bind (a b)
                (skip-loop num skip-json loop-json staple-vec scaffold-vec)
              (setf (p-strand vstrand) a
                    (n-strand vstrand) b)))
      ;;fill base names
      (format t "About to apply node-staple-sequences~%")
      (loop for sequence being the hash-values in node-staple-sequences using (hash-key start-node)
         do (loop with node = start-node
               for base across sequence
               do (unless (char= base #\?)
                    (setf (base-name node) (intern base :keyword)))
                 (setf node (forward-node node))
                 (unless node (return))))
      vstrands)))

   
(defun read-csv-sequence-file (csv-filename)
  (if (probe-file csv-filename)
      (with-open-file (csv-filename :direction :input)
        (read-csv:parse-csv csv-file))
      nil))

(defun find-title-row (data)
  (find-if (lambda (row)
             (alpha-char-p (char (first row) 0)))
           data))

(defun title-positions (title-row)
  (values (position "Start" title-row :test #'string=)
          (position "Sequence" title-row :test #'string=)))

(defun positions-and-sequences (csv)
  (let* ((title-row (find-title-row csv))
         (data (remove title-row csv)))
    (multiple-value-bind (start sequence)
        (title-positions title-row)
      (loop for row in data
         collect (multiple-value-bind (strand pos)
                     (position-from-string (nth start row))
                   (list strand pos
                         (nth sequence row)))))))

(defun position-from-string (string)
  (multiple-value-bind (first stop)
      (parse-integer string :junk-allowed t)
    (values first (parse-integer string :start (1+ stop) :junk-allowed t))))

(defun apply-default-sequence-to-cstrands (cstrands)
  (loop for strand in cstrands
     do (loop for node in (dna-strand-as-list-of-nodes strand)
           do (progn
                (setf (base-name node) :G)
                (when (hbond-node node)
                  (setf (base-name (hbond-node node)) :C))))))

(defun node-has-identifier-p (node)
  (slot-boundp node 'identifier))

(defun node-has-base-name-p (node)
  (slot-boundp node 'base-name))

(defun complementary-base (base)
  (ecase base
    ((:c) :g)
    ((:g) :c)
    ((:a) :t)
    ((:t) :a)))

(defun fill-empty-bases-with-default (cstrands)
  ;;; Apply the sequence information in csv-sequence to the vstrands
  ;;;    If there is no csv-sequence then use :G/:C
  (loop for strand in cstrands
     do (loop for node in (dna-strand-as-list-of-nodes strand)
           for hnode = (hbond-node node)
           do (format t "fill-empty-bases-with-default node -> ~a~%" node)
           do (if (node-has-base-name-p node)
                  (progn
                    (format t "TRUE node-has-base-name-p~%")
                    (cond
                      ((not hnode))
                      ((node-has-base-name-p hnode)
                       (unless (eql (base-name hnode)
                                    (complementary-base (base-name node)))
                         (error "Base pair doesn't match: ~a" node)))
                      (t (setf (base-name hnode)
                               (complementary-base (base-name node))))))
                  (progn
                    (format t "NOT node-has-base-name-p~%")
                    (cond ((not hnode)
                         (setf (base-name node) :G))
                        ((node-has-base-name-p hnode)
                         (setf (base-name node)
                               (complementary-base (base-name hnode))))
                        (t (setf (base-name node) :G
                                 (base-name hnode) :C))))))))

(defun parse-cadnano (file-name)
  (let ((json-filename (make-pathname :type "json" :defaults file-name))
        (csv-filename (make-pathname :type "csv" :defaults file-name)))
    (let* ((json (with-open-file (json-file json-filename :direction :input)
                   (json:decode-json json-file)))
           (csv (with-open-file (csv-file csv-filename :direction :input
                                          :if-does-not-exist nil)
                  (positions-and-sequences
                   (read-csv:parse-csv csv-file))))
           (vstrands (parse-json json csv))
           (cstrands (single-or-double-classified-strands vstrands)))
      (format t "About to fill-empty-bases-with-default~%")
      (fill-empty-bases-with-default cstrands)
      cstrands)))
    
(defun BUILD-NODE (strand-json vec num strand-name)
  (loop for index from 0 below (length vec)
     for node-json in strand-json
     unless (apply #'= -1 node-json)
     do (setf (elt vec index) (make-instance 'node :identifier (list :j num strand-name index)))))

(defun INTRA-HELIX-CONNECT (staple-vec scaffold-vec)
  (loop for index from 0 below (length staple-vec)
     for staple-node = (elt staple-vec index)
     for scaffold-node = (elt scaffold-vec index)
     do (when (and staple-node scaffold-node)
          (setf (hbond-node staple-node) scaffold-node)
          (setf (hbond-node scaffold-node) staple-node))))
    
(defun connect-everything (vstrands num strand-json strand-vec accessor)
  (loop for json-info in strand-json
     for index from 0 below (length strand-json)
     for node = (elt strand-vec index)
     when node
     do (let* ((fwd-vstrand (first json-info))
               (fwd-pos     (second json-info))
               (bkd-vstrand (third json-info))
               (bkd-pos     (fourth json-info)))
          (when (/= fwd-vstrand -1)
            (let ((forward-node (lookup-node vstrands fwd-vstrand fwd-pos accessor)))
              (assert forward-node)
              (setf (forward-node node) forward-node)))
          (when (/= bkd-vstrand -1)
            (let ((backward-node (lookup-node vstrands bkd-vstrand bkd-pos accessor)))
              (assert backward-node)
              (setf (backward-node node) backward-node))))))
              
(defclass node ()
  ((hbond-node :initform nil :initarg :hbond-node :accessor hbond-node)
   (t-q-bond-node :initform nil :initarg :t-q-bond-node :accessor t-q-bond-node)
   (forward-node :initform nil :initarg :forward-node :accessor forward-node)
   (backward-node :initform nil :initarg :backward-node :accessor backward-node)
   (base-name :initarg :base-name :accessor base-name)
   (residue :initform nil :initarg :residue :accessor residue)
   (strand :initarg :strand :accessor strand)
   (identifier :initarg :identifier :reader identifier)))

(defun lookup-node (vstrands vstrand-num pos accessor)
  (let* ((vstrand (gethash vstrand-num vstrands))
         (strand (funcall accessor vstrand)))
    (elt strand pos)))


(defun arrow-direction (vec)
  (let ((step-direction (loop for index from 0 below (length vec)
                           for node = (elt vec index)
                           when node
                           collect (let ((fwd (forward-node node)))
                                     (when (and fwd (position fwd vec))
                                       (- (position fwd vec) index))))))
#|    (cond
      ((every (lambda (x) (if x (plusp x) T)) step-direction) 1)
      ((every (lambda (x) (if x (minusp x) T)) step-direction) -1)
      (t (error "arrow direction in the vector is neither forward nor backward ~a" step-direction)))))|#
                                        ;(if (every nil step-direction) (error "empty vector")
    (if step-direction 
        (loop for step in step-direction
           with pos = 0 and neg = 0
           when step
           do (cond ((plusp step) (incf pos))
                    ((minusp step) (incf neg))
                    (T (error "Do we find a basis connected with itself? ~a" step-direction)))
           finally (return (cond ((> pos neg) 1)
                                 ((< pos neg) -1)
                                 (T (error "arrow direction in the vector is neither forward nor backward ~a" step-direction)))))
        nil)))

(defun skip-loop (num skip-json loop-json staple-vec scaffold-vec)
  (flet ((skip-procedure (x)
           (let ((xf (forward-node x))
                 (xb (backward-node x)))
             (when xb (setf (forward-node xb) xf))
             (when xf (setf (backward-node xf) xb))
             (setf x NIL)))
         (loop-procedure (num strand-name old-vec s-l-cursor new-vec dest-cursor loop-vec setf-accessor-l accessor-l setf-accessor-r)
            (setf (elt new-vec dest-cursor) (make-instance 'node :identifier (list :A num strand-name s-l-cursor :I (elt loop-vec s-l-cursor))))
            (let ((New-node (elt new-vec dest-cursor))
                  (Doubled-node (elt old-vec s-l-cursor))
                  (l-d-node (funcall accessor-l (elt old-vec s-l-cursor))))
              (when Doubled-node (funcall setf-accessor-l New-node Doubled-node))
              (funcall setf-accessor-r  Doubled-node New-node)
              (funcall setf-accessor-l l-d-node New-node)
              (when l-d-node (funcall setf-accessor-r New-node l-d-node)))))
    ;;must check if 4 input have the same lenght
    (let ((skip-vec (coerce skip-json 'vector))
          (loop-vec (coerce loop-json 'vector))
          (stap-direction (arrow-direction staple-vec))
          (scaf-direction (arrow-direction scaffold-vec)))
      (format t "Starting function~%")
      (let* ((old-vec-min-length
              (loop for index from 0 below (length scaffold-vec)
                 when (or (elt scaffold-vec index) (elt staple-vec index))
                 count 1))
             (new-vec-length (+ old-vec-min-length (reduce #'+ skip-json) (reduce #'+ loop-json)))
             (p-strand (make-array new-vec-length));p-strand
             (n-strand (make-array new-vec-length))
             (condition (cond
                          ((and (not stap-direction)
                                (not scaf-direction))
                           t)
                          ((and (or (not stap-direction) (= stap-direction 1))
                                (or (not scaf-direction) (= scaf-direction -1)))
                           t)
                          ((and (or (not stap-direction) (= stap-direction -1))
                                (or (not scaf-direction) (= scaf-direction 1)))
                           nil)
                          (t (error "What do I do with stap-direction = ~a and scaf-direction = ~a~%" stap-direction scaf-direction))))
             (old-p-strand (if condition staple-vec scaffold-vec))
             (old-n-strand (if condition scaffold-vec staple-vec)))
        (format t "Starting loop~%")
        (loop with s-l-cursor = 0
           with dest-cursor = 0
           with length-skip = (length skip-vec)
           until (= s-l-cursor length-skip)
           do (progn
;               (format t "s-l-cursor -> ~a~%" s-l-cursor)
                (if (or (elt scaffold-vec s-l-cursor) (elt staple-vec s-l-cursor))
                    (cond
                      ((and (= (elt skip-vec s-l-cursor) -1)
                            (= (elt loop-vec s-l-cursor) 0))
                       ;;do a skip procedure
                       (format t "hit a skip~%")
                       (skip-procedure (elt staple-vec s-l-cursor))
                       (skip-procedure (elt scaffold-vec s-l-cursor))
                       (incf s-l-cursor)
                       )
                      ((and (= (elt skip-vec s-l-cursor) 0)
                            (> (elt loop-vec s-l-cursor) 0))
                       ;;do an insertion
                       (format t "Insertion number -> ~a~%" (elt loop-vec s-l-cursor))
                       (and (elt staple-vec s-l-cursor)
                            (loop-procedure num :ST old-p-strand s-l-cursor p-strand dest-cursor loop-vec #'(setf backward-node) #'backward-node #'(setf forward-node))); p-strand
                       (and (elt scaffold-vec s-l-cursor)
                            (loop-procedure num :SC old-n-strand s-l-cursor n-strand dest-cursor loop-vec #'(setf forward-node) #'forward-node #'(setf backward-node))); n-strand
                      (let ((p-node (elt p-strand dest-cursor)); p-node p-strand
                             (n-node (elt n-strand dest-cursor))); n-node n-strand
                         (when (and p-node n-node); p-node n-node
                           (setf (hbond-node p-node) n-node); p-node n-node
                           (setf (hbond-node n-node) p-node))); n-node p-node 
                       (decf (elt loop-vec s-l-cursor))
                       (incf dest-cursor))
                      ((and (= (elt skip-vec s-l-cursor) 0)
                            (= (elt loop-vec s-l-cursor) 0))
                       ;;do a copy
;                      (format t "copy~%")
                       (setf (elt p-strand dest-cursor) (elt old-p-strand s-l-cursor)); n or p
                       (setf (elt n-strand dest-cursor) (elt old-n-strand s-l-cursor));n or p
                       (incf dest-cursor)
                       (incf s-l-cursor))
                      (t (error "An impossible step was encountered")))
                    (incf s-l-cursor)
                    )))
        (values p-strand n-strand); n-strand p-strand
;       (list new-stap-vec new-scaf-vec)
;       new-stap-vec
        ))))


