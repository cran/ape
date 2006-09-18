### all.equal.phylo.R  (2006-09-12)
###
###     Global Comparison of two Phylogenies
###
### Copyright 2006 Benoît Durand
###
### This file is part of the `ape' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

### Recherche de la correspondance entre deux arbres
### Parcours en profondeur et en parallèle des deux arbres (current et target)
### current, target: les deux arbres à comparer
### use.edge.length: faut-il comparer les longueurs de branches ?
### use.tip.label: faut-il comparer les étiquettes de feuilles ou seulement la
###	topologie des deux arbres ?
### index.return: si TRUE, retourner la matrice de correspondance entre noeuds
###	et feuilles, une matrice à deux colonnes (current et target) avec pour
###	chaque ligne des paires d'identifiants de noeuds/feuilles, tels qu'ils
###	apparaissent dans l'attribut 'edge' des objets phylo
### tolerance, scale: paramètres de comparaison des longueurs de branches
###	(voir 'all.equal')
all.equal.phylo <- function(target, current,
                        use.edge.length = TRUE,
                        use.tip.label = TRUE,
                        index.return = FALSE,
                        tolerance = .Machine$double.eps ^ 0.5,
                        scale = NULL, ...) {

	same.node <- function(i, j) {
		ni <- as.numeric(i)
		nj <- as.numeric(j)
		# Comparaison de un noeud et une feuille
		if (sign(ni) != sign(nj)) return(NULL)
		# Comparaison de deux feuilles
		if (sign(ni) == 1) {
			if (!use.tip.label) return(c(i, j))
			if (current$tip.label[ni] == target$tip.label[nj])
				return(c(i, j))
			return(NULL)
  		}
  		# Comparaison de deux noeuds
		i.children <- which(current$edge[, 1] == i)
		j.children <- which(target$edge[, 1] == j)
		if (length(i.children) != length(j.children)) return(NULL)
		correspondance <- NULL
		for (i.child in i.children) {
			corresp <- NULL
			for (j.child in j.children) {
				if (!use.edge.length ||
					isTRUE(all.equal(current$edge.length[i.child],
							     target$edge.length[j.child],
							     tolerance = tolerance,
							     scale = scale)))
					corresp <- same.node(current$edge[i.child, 2],
								   target$edge[j.child, 2])
				if (!is.null(corresp)) break
			}
			if (is.null(corresp)) return(NULL)
			correspondance <- c(correspondance, i, j, corresp)
			j.children <- j.children[j.children != j.child]
		}
		return(correspondance)
	}

	result <- same.node('-1','-1')
	if (!isTRUE(index.return)) return(!is.null(result))
	if (is.null(result)) return(result)
	result <- t(matrix(result, nrow = 2))
      colnames(result) = c('current', 'target')
      return(result)
}
