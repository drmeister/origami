(in-package :origami)

;;; ------------------------------------------------------------
;;;
;;;  Graphviz generated

(defun safe-draw-link (dest node next-node color &optional link)
  ;(let ((*print-pretty* nil)))
  (when (and node next-node)
    (format dest "\"~a\" -> \"~a\" [color=~a];~%" (identifier node) (identifier next-node) color)))

(defun draw-connection (dest node)
  (safe-draw-link dest node (forward-node node) "violet")
  (safe-draw-link dest node (backward-node node) "blue")
  (safe-draw-link dest node (hbond-node node) "orange"))

(defun draw-node (dest node)
;  (let ((*print-pretty* nil)))
  (format dest "\"(~{~a~^ ~})\"->" (identifier node)))

(defun draw-strand (dest num name vec)
  (format dest "subgraph cluster~a_~a {~%" (string name) num )
  (format dest "     edge [style=\"invis\"];~%")
  (format dest "     label = \"~a\";~%" (string name))
  (format dest "     color = blue;~%")
  (format dest "     \"start~a_~a\"  [style=\"invis\"];~%" (string name) num)
  (format dest "     \"stop~a_~a\"  [style=\"invis\"];~%" (string name) num)
  (format dest "     \"start~a_~a\"->" (string name) num)
  ;(direction (arrow-direction vec))
  (loop for node across vec
     ;for direction =(arrow-direction vec); (if (= 1 (arrow-direction vec)) (#'forward-node) (#'backward-node))
     when node
     do (draw-node dest node))
  (format dest "\"stop~a_~a\"}~%" (string name) num))
  
(defun draw-double-strand (dest num scaffold staple)
  (format dest "subgraph cluster_~a {~%" num)
  (format dest "label = \"double-strand-~a\";~%" num)
  (draw-strand dest num :scaffold scaffold)
  (draw-strand dest num :staple staple)
  (format dest "}~%"))


(defun draw-graph (dest strands)
  (let ((dest (if (streamp dest)
                  dest
                  (open dest :direction :output :if-exists :supersede))))
    (format dest "digraph G {~%")
    (loop for strand being the hash-values in strands using (hash-key num)
       for staple-vec = (staple-vec strand)
       for scaffold-vec = (scaffold-vec strand)
       do (draw-double-strand dest num scaffold-vec staple-vec))
    (loop for strand being the hash-values in strands using (hash-key num)
       for staple-vec = (staple-vec strand)
       for scaffold-vec = (scaffold-vec strand)
       do (loop for node across scaffold-vec
             when node
             do (draw-connection dest node))
       do (loop for node across staple-vec
             when node
             do (draw-connection dest node))
         )
    (format dest "}~%")))

(defun draw-result-double-strand (dest num p-strand n-strand)
  (format dest "subgraph cluster_~a {~%" num)
  (format dest "label = \"double-strand-~a\";~%" num)
  (draw-strand dest num :p p-strand)
  (draw-strand dest num :n n-strand)
  (format dest "}~%"))

(defun draw-result-graph (dest strands)
  (let ((dest (if (streamp dest)
                  dest
                  (open dest :direction :output :if-exists :supersede))))
;    (format t "1")
    (format dest "digraph G {~%")
    (loop for strand being the hash-values in strands using (hash-key num)
       for p-strand = (p-strand strand)
       for n-strand = (n-strand strand)
;       do (format t "2")
       do (draw-double-strand dest num p-strand n-strand))
    (loop for strand being the hash-values in strands using (hash-key num)
       for p-strand = (p-strand strand)
       for n-strand = (n-strand strand)
       do (loop for node across p-strand
             when node
             do (draw-connection dest node))
       do (loop for node across n-strand
             when node
             do (draw-connection dest node))
         )
    (format dest "}~%")))

