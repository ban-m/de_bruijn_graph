mod find_union;
use find_union::FindUnion;
use log::debug;
use std::collections::{HashMap, HashSet};
#[derive(Clone)]
pub struct DeBruijnGraph {
    pub k: usize,
    pub nodes: Vec<Node>,
    indexer: HashMap<Node, usize>,
}

impl std::fmt::Debug for DeBruijnGraph {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "NumOfNodes:{}", self.nodes.len())?;
        for (idx, node) in self.nodes.iter().enumerate() {
            writeln!(f, "{}\t{:?}", idx, node)?;
        }
        write!(f, "K:{}", self.k)
    }
}

#[derive(Clone)]
pub struct Node {
    pub occ: usize,
    pub edges: Vec<Edge>,
    pub kmer: Vec<(u64, u64)>,
    pub cluster: Option<usize>,
}

impl Node {
    fn push(&mut self, to: usize) {
        match self.edges.iter_mut().find(|e| e.to == to) {
            Some(x) => {
                x.weight += 1;
            }
            None => self.edges.push(Edge { to, weight: 1 }),
        }
    }
    fn remove_edge(&mut self, to: usize) {
        self.edges.retain(|x| x.to != to);
    }
    pub fn new(kmer: Vec<(u64, u64)>) -> Self {
        Node {
            kmer: kmer.to_vec(),
            edges: vec![],
            occ: 0,
            cluster: None,
        }
    }
}

impl std::fmt::Debug for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let edges: Vec<_> = self
            .edges
            .iter()
            .map(|e| format!("(->{},{})", e.to, e.weight))
            .collect();
        let kmer: Vec<_> = self
            .kmer
            .iter()
            .map(|(u, c)| format!("{}-{}", u, c))
            .collect();
        let cluster = self
            .cluster
            .map(|x| format!("{}", x))
            .unwrap_or_else(|| "-".to_string());
        write!(
            f,
            "{}\t{}\t{}\t[{}]",
            cluster,
            self.occ,
            edges.join(","),
            kmer.join(",")
        )
    }
}

pub trait AsDeBruijnNode: std::marker::Sized {
    fn as_node(w: &[Self]) -> Node;
}

#[derive(Debug, Clone)]
pub struct Edge {
    to: usize,
    pub weight: u64,
}

use std::hash::Hasher;
impl std::hash::Hash for Node {
    fn hash<H: Hasher>(&self, state: &mut H) {
        assert!(!self.kmer.is_empty());
        // Normalize and hashing.
        if self.kmer.first().unwrap() < self.kmer.last().unwrap() {
            for (unit, cluster) in self.kmer.iter() {
                unit.hash(state);
                cluster.hash(state);
            }
        } else {
            for (unit, cluster) in self.kmer.iter().rev() {
                unit.hash(state);
                cluster.hash(state);
            }
        }
    }
}

impl std::cmp::PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        assert!(!self.kmer.is_empty());
        assert!(!other.kmer.is_empty());
        if self.kmer.len() != other.kmer.len() {
            return false;
        }
        let is_self_normed = self.kmer.first().unwrap() < self.kmer.last().unwrap();
        let is_other_normed = other.kmer.first().unwrap() < other.kmer.last().unwrap();
        match (is_self_normed, is_other_normed) {
            (false, false) | (true, true) => self.kmer.iter().zip(&other.kmer).all(|(x, y)| x == y),
            (false, true) | (true, false) => {
                self.kmer.iter().rev().zip(&other.kmer).all(|(x, y)| x == y)
            }
        }
    }
}
impl std::cmp::Eq for Node {}

pub trait IntoDeBruijnNodes {
    fn into_de_bruijn_nodes(&self, k: usize) -> Vec<Node>;
}

impl DeBruijnGraph {
    pub fn from<T: IntoDeBruijnNodes>(reads: &[T], k: usize) -> Self {
        let (mut nodes, mut indexer) = (vec![], HashMap::new());
        for read in reads {
            let read = read.into_de_bruijn_nodes(k);
            for (idx, w) in read.windows(2).enumerate() {
                let (from, to) = (&w[0], &w[1]);
                // Check entry.
                let from = if !indexer.contains_key(from) {
                    indexer.insert(from.clone(), nodes.len());
                    nodes.push(from.clone());
                    nodes.len() - 1
                } else {
                    *indexer.get(from).unwrap()
                };
                let to = if !indexer.contains_key(to) {
                    indexer.insert(to.clone(), nodes.len());
                    nodes.push(to.clone());
                    nodes.len() - 1
                } else {
                    *indexer.get(to).unwrap()
                };
                nodes[from].push(to);
                nodes[to].push(from);
                if idx == 0 {
                    nodes[from].occ += 1;
                }
                nodes[to].occ += 1;
            }
        }
        Self { k, nodes, indexer }
    }
    pub fn validate(&self) {
        for idx in 0..self.nodes.len() {
            for edge in self.nodes[idx].edges.iter() {
                assert!(self.nodes[edge.to].edges.iter().any(|e| e.to == idx));
            }
        }
    }
    #[allow(dead_code)]
    fn calc_thr_edge(&self) -> u64 {
        let counts = self
            .nodes
            .iter()
            .map(|n| n.edges.iter().fold(0, |x, w| x + w.weight))
            .sum::<u64>();
        let len: usize = self.nodes.iter().map(|n| n.edges.len()).sum::<usize>();
        counts / len as u64 / 3
    }
    #[allow(dead_code)]
    fn calc_thr(&self) -> u64 {
        let mut counts: HashMap<Vec<u64>, usize> = HashMap::new();
        for node in self.nodes.iter() {
            let kmer: Vec<_> = node.kmer.iter().map(|n| n.0).collect();
            *counts.entry(kmer).or_default() += node.occ;
        }
        let mut hist: HashMap<_, u32> = HashMap::new();
        for count in counts.values() {
            *hist.entry(*count).or_default() += 1;
        }
        let mut probe = 1;
        while hist.contains_key(&probe) {
            probe += 1;
        }
        // Changed.
        3 * probe as u64 - 1
    }
    pub fn remove_isolated_nodes(&mut self) {
        let to_be_removed: Vec<_> = self.nodes.iter().map(|e| e.edges.is_empty()).collect();
        self.remove_nodes(&to_be_removed);
    }
    pub fn clean_up_tight(&mut self) {
        let total = self.nodes.iter().map(|n| n.edges.len()).sum::<usize>();
        let len = self.nodes.len();
        for i in 0..len {
            if self.nodes[i].edges.len() <= 2 {
                continue;
            }
            self.clean_up_tight_at(i);
        }
        let after = self.nodes.iter().map(|n| n.edges.len()).sum::<usize>();
        debug!("{}=>{}", total, after);
    }
    fn clean_up_tight_at(&mut self, i: usize) {
        // Enumerate candidates.
        let mut degree: HashMap<u64, u32> = HashMap::new();
        let kmer = self.nodes[i].kmer.clone();
        for edge in self.nodes[i].edges.iter() {
            for &(pos, cl) in self.nodes[edge.to].kmer.iter() {
                if !kmer.contains(&(pos, cl)) {
                    *degree.entry(pos).or_default() += 1;
                }
            }
        }
        for (pos, _) in degree.into_iter().filter(|&(_, c)| c > 1) {
            // Select edges towards pos.
            let max = self.nodes[i]
                .edges
                .iter()
                .filter(|e| self.nodes[e.to].kmer.iter().any(|x| x.0 == pos))
                .map(|e| e.weight)
                .max()
                .unwrap();
            let removed: Vec<usize> = self.nodes[i]
                .edges
                .iter()
                .filter(|w| {
                    let has_pos = self.nodes[w.to].kmer.iter().any(|x| x.0 == pos);
                    let heavy = w.weight > max / 2;
                    has_pos && !heavy
                })
                .map(|e| e.to)
                .collect();
            for r in removed {
                self.nodes[i].remove_edge(r);
                self.nodes[r].remove_edge(i);
            }
        }
    }
    pub fn clean_up_auto(&mut self) {
        let thr = self.calc_thr();
        debug!("Removing edges with weight less than {}", thr);
        self.clean_up(thr)
    }
    pub fn clean_up(&mut self, thr: u64) {
        let mut removed_edge = vec![];
        self.nodes
            .iter_mut()
            .enumerate()
            .filter(|(_, n)| n.edges.len() > 2)
            .for_each(|(idx, node)| {
                node.edges.retain(|edge| {
                    if edge.weight > thr {
                        true
                    } else {
                        removed_edge.push((idx, edge.to));
                        false
                    }
                })
            });
        for (to, from) in removed_edge {
            self.nodes[from].remove_edge(to);
        }
        self.validate();
    }
    pub fn assign_read_by_unit<T: IntoDeBruijnNodes>(&self, read: &T) -> Option<usize> {
        let max_cluster = self.nodes.iter().flat_map(|n| n.cluster).max().unwrap_or(0);
        let mut units_in_cluster: Vec<HashSet<_>> = vec![HashSet::new(); max_cluster + 1];
        for node in self.nodes.iter() {
            if let Some(cl) = node.cluster {
                for tuple in node.kmer.iter() {
                    units_in_cluster[cl as usize].insert(tuple.clone());
                }
            }
        }
        units_in_cluster
            .iter()
            .enumerate()
            .map(|(idx, set)| {
                let count = read
                    .into_de_bruijn_nodes(1)
                    .iter()
                    .filter(|node| set.contains(&node.kmer[0]))
                    .count();
                (idx, count)
            })
            .filter(|x| x.1 > 0)
            .max_by_key(|x| x.1)
            .map(|x| x.0)
    }
    pub fn assign_read<T: IntoDeBruijnNodes>(&self, read: &T) -> Option<usize> {
        let mut count = HashMap::<_, u32>::new();
        for node in read.into_de_bruijn_nodes(self.k) {
            if let Some(&idx) = self.indexer.get(&node) {
                if let Some(cl) = self.nodes[idx].cluster {
                    *count.entry(cl).or_default() += 1;
                }
            }
        }
        count.into_iter().max_by_key(|x| x.1).map(|x| x.0)
    }
    pub fn coloring_nodes_by<T: IntoDeBruijnNodes>(&mut self, read: &T, cluster: usize) {
        for node in read.into_de_bruijn_nodes(self.k) {
            if let Some(&idx) = self.indexer.get(&node) {
                self.nodes[idx].cluster = Some(cluster);
            }
        }
    }
    pub fn expand_color<T: IntoDeBruijnNodes>(
        &mut self,
        reads: &[T],
        thr: usize,
        labels: &[Option<u8>],
        forbs: &[&[u8]],
        background: Option<u8>,
    ) -> HashMap<u8, u8> {
        debug!("{:?}", background);
        let mut current_cluster =
            labels.iter().cloned().filter_map(|x| x).max().unwrap_or(0) as usize + 1;
        let mut fu = FindUnion::new(self.nodes.len());
        for (from, node) in self.nodes.iter().enumerate().filter(|x| x.1.occ > 0) {
            for edge in node.edges.iter().filter(|e| e.weight > 0) {
                fu.unite(from, edge.to);
            }
        }
        let mut map = HashMap::new();
        debug!("ClusterID\tNumberOfKmer");
        let mut ignored = 0;
        for cluster in 0..self.nodes.len() {
            if fu.find(cluster).unwrap() != cluster {
                continue;
            }
            let count = (0..self.nodes.len())
                .filter(|&x| fu.find(x).unwrap() == cluster)
                .count();
            if count < thr {
                ignored += 1;
                continue;
            }
            let mut counts: HashMap<_, u32> = HashMap::new();
            let mut forbidden_cluster: HashMap<_, u32> = HashMap::new();
            for (i, read) in reads.iter().enumerate() {
                let (total, contained) =
                    read.into_de_bruijn_nodes(self.k)
                        .iter()
                        .fold((0, 0), |(total, cont), node| {
                            let is_in_cluster = match self.indexer.get(&node) {
                                Some(&idx) => (fu.find(idx).unwrap() == cluster) as u32,
                                None => 0,
                            };
                            (total + 1, cont + is_in_cluster)
                        });
                if total / 3 < contained {
                    if let Some(label) = labels[i] {
                        *counts.entry(label).or_default() += 1;
                    }
                    for f in forbs[i].iter() {
                        *forbidden_cluster.entry(f).or_default() += 1;
                    }
                }
            }
            let merged_clusters: Vec<_> = counts
                .iter()
                .filter(|&(_, &count)| count > 10)
                .map(|c| c.0)
                .collect();
            let assignment = match merged_clusters.iter().min() {
                Some(&&res) => res as usize,
                None if background.is_some() => background.unwrap() as usize,
                None => {
                    current_cluster += 1;
                    current_cluster - 1
                }
            };
            for (idx, node) in self.nodes.iter_mut().enumerate() {
                if fu.find(idx).unwrap() == cluster {
                    node.cluster = Some(assignment);
                }
            }
            for &&x in merged_clusters.iter() {
                map.insert(x, assignment as u8);
            }
            debug!(
                "{}\t{}\t{:?}\t{:?}",
                assignment, count, forbidden_cluster, merged_clusters
            );
        }
        debug!("Ignored small component (<{} kmer):{}", thr, ignored);
        {
            let mut res: HashMap<_, u32> = HashMap::new();
            for node in self.nodes.iter() {
                if let Some(cl) = node.cluster {
                    *res.entry(cl).or_default() += 1;
                }
            }
            debug!("Cluster:{:?}", res);
        }
        map
    }
    pub fn coloring(&mut self, thr: usize) -> usize {
        // Coloring node of the de Bruijn graph.
        // As a first try, I just color de Bruijn graph by its connected components.
        let mut fu = FindUnion::new(self.nodes.len());
        for (from, node) in self.nodes.iter().enumerate().filter(|x| x.1.occ > 0) {
            for edge in node.edges.iter().filter(|e| e.weight > 0) {
                fu.unite(from, edge.to);
            }
        }
        let mut current_component = 0;
        debug!("ClusterID\tNumberOfKmer");
        let mut ignored = 0;
        for cluster in 0..self.nodes.len() {
            if fu.find(cluster).unwrap() != cluster {
                continue;
            }
            let count = (0..self.nodes.len())
                .filter(|&x| fu.find(x).unwrap() == cluster)
                .count();
            if count < thr {
                ignored += 1;
                continue;
            }
            debug!("{}\t{}", current_component, count);
            for (idx, node) in self.nodes.iter_mut().enumerate() {
                if fu.find(idx).unwrap() == cluster {
                    node.cluster = Some(current_component);
                }
            }
            current_component += 1;
        }
        debug!("Ignored small component (<{} kmer):{}", thr, ignored);
        {
            let mut res: HashMap<_, u32> = HashMap::new();
            for node in self.nodes.iter() {
                if let Some(cl) = node.cluster {
                    *res.entry(cl).or_default() += 1;
                }
            }
            debug!("Cluster:{:?}", res);
        }
        current_component
    }
    #[allow(dead_code)]
    pub fn clustering(&self, thr: usize) -> Vec<HashSet<(u64, u64)>> {
        // Clustering de Bruijn graph.
        // As a first try, I implement very naive conneceted component analysis.
        // To this end, I use naive FindUnion Tree. In other words,
        // as I traverse nodes, I merge the two connected nodes.
        // Currently, I ignore very weak connection, i.e., connection
        // with the weight of less than 1.
        let mut fu = FindUnion::new(self.nodes.len());
        for (from, node) in self.nodes.iter().enumerate().filter(|(_, n)| n.occ > thr) {
            for edge in node.edges.iter().filter(|e| e.weight > thr as u64) {
                fu.unite(from, edge.to);
            }
        }
        let mut components = vec![];
        for cluster in 0..self.nodes.len() {
            if fu.find(cluster).unwrap() != cluster {
                continue;
            }
            let count = (0..self.nodes.len())
                .filter(|&x| fu.find(x).unwrap() == cluster)
                .count();
            debug!("Find cluster. Containing {} k-mers.", count);
            let component: HashSet<_> = (0..self.nodes.len())
                .filter(|&x| fu.find(x).unwrap() == cluster)
                .flat_map(|node_idx| self.nodes[node_idx].kmer.iter())
                .copied()
                .collect();
            if component.len() > 1 {
                components.push(component);
            }
        }
        components
    }
    /// Resolve bubbles in this graph. Note that after calling this function,
    /// Some internal nodes would be removed, and thus the connection
    /// between two k-mers might be invalid. Specifically, all node between
    /// two bubbles would be removed if the bubble can be merged.
    pub fn resolve_bubbles<T: IntoDeBruijnNodes>(&mut self, reads: &[T]) {
        let mut bubble_spec = self.enumerate_bubbles(reads);
        let mut queue_and_parent: std::collections::VecDeque<_> = bubble_spec
            .iter()
            .enumerate()
            .filter_map(|(idx, e)| e.map(|_| (idx, idx)))
            .collect();
        let mut to_be_removed = vec![false; self.nodes.len()];
        while let Some((idx, parent)) = queue_and_parent.pop_back() {
            let bubble = match bubble_spec[idx] {
                Some(x) => x,
                None => continue,
            };
            let pair_pos = match self.search_nearest_bubble(bubble.root, bubble.shoot, &bubble_spec)
            {
                Ok(res) => res,
                Err(ReachedNodeType::Root) => {
                    // We cut the bubble arbitrary.
                    let branch = bubble.branches.0;
                    self.nodes[branch].remove_edge(bubble.shoot);
                    self.nodes[bubble.shoot].remove_edge(branch);
                    continue;
                }
                Err(ReachedNodeType::BubbleBranch(p)) if parent == p => {
                    // We reached the same parent again.
                    // I think this bubble could not be resolved. Just abandon.
                    continue;
                }
                Err(ReachedNodeType::BubbleBranch(p)) => {
                    // We reached a branch of a bubble. However, thre is
                    // some progress.
                    // Hope we can resolve this bubble in the next time.
                    queue_and_parent.push_front((idx, p));
                    continue;
                }
                Err(ReachedNodeType::ComplexNode) => continue,
            };
            let pair = bubble_spec[pair_pos].unwrap();
            bubble_spec[idx] = None;
            bubble_spec[pair_pos] = None;
            // contract into bubble.shoot.
            {
                let mut child = vec![];
                self.nodes[pair.branches.0]
                    .edges
                    .iter_mut()
                    .filter(|e| e.to == pair.shoot)
                    .for_each(|e| {
                        child.push((pair.branches.0, e.weight));
                        e.to = bubble.shoot
                    });
                self.nodes[pair.branches.1]
                    .edges
                    .iter_mut()
                    .filter(|e| e.to == pair.shoot)
                    .for_each(|e| {
                        child.push((pair.branches.1, e.weight));
                        e.to = bubble.shoot;
                    });
                for (to, w) in child {
                    let edges = &mut self.nodes[bubble.shoot].edges;
                    match edges.iter_mut().find(|e| e.to == to) {
                        Some(x) => {
                            x.weight += w;
                        }
                        None => edges.push(Edge { to, weight: w }),
                    }
                }
            }
            // Register all the nodes between two shoots as `removed`
            let (mut prev, mut current) = (bubble.shoot, bubble.root);
            while current != pair.shoot {
                to_be_removed[current] = true;
                let next_nodes = &self.nodes[current].edges;
                assert!(next_nodes.len() == 2);
                let next = next_nodes.iter().find(|e| e.to != prev).unwrap().to;
                self.nodes[current].edges.clear();
                prev = current;
                current = next;
            }
            self.nodes[pair.shoot].edges.clear();
            self.nodes[bubble.shoot].remove_edge(bubble.root);
            to_be_removed[pair.shoot] = true;
            // solve crossing.
            let merged_nodes = match self.resolve_crossing(bubble.shoot, reads) {
                Some(res) => res,
                None => continue,
            };
            // Success.
            to_be_removed[bubble.shoot] = true;
            let children: Vec<_> = self.nodes[bubble.shoot]
                .edges
                .iter()
                .map(|e| e.to)
                .collect();
            for child in children {
                self.nodes[child].remove_edge(bubble.shoot);
            }
            self.nodes[bubble.shoot].edges.clear();
            for (contract, members) in merged_nodes {
                for member in members {
                    self.nodes[member].push(contract);
                    self.nodes[contract].push(member);
                }
            }
        }
        // Removing nodes.
        self.remove_nodes(&to_be_removed);
    }
    fn enumerate_bubbles<T: IntoDeBruijnNodes>(&self, reads: &[T]) -> Vec<Option<Bubble>> {
        // Enumerate bubble.
        let mut edge_counts: Vec<_> = (0..self.nodes.len()).map(|idx| vec![0; idx]).collect();
        for read in reads.iter() {
            let read = read.into_de_bruijn_nodes(self.k);
            for w in read.windows(3) {
                let from = match self.indexer.get(&w[0]) {
                    Some(&res) => res,
                    None => continue,
                };
                let to = match self.indexer.get(&w[2]) {
                    Some(&res) => res,
                    None => continue,
                };
                edge_counts[from.max(to)][from.min(to)] += 1;
            }
        }
        self.nodes
            .iter()
            .enumerate()
            .map(|(shoot, node)| {
                if node.edges.len() != 3 {
                    return None;
                }
                let (to0, to1, to2) = (node.edges[0].to, node.edges[1].to, node.edges[2].to);
                let bet0and1 = edge_counts[to0.max(to1)][to0.min(to1)];
                let bet0and2 = edge_counts[to0.max(to2)][to0.min(to2)];
                let bet1and2 = edge_counts[to1.max(to2)][to1.min(to2)];
                let (branches, root) = if bet0and1 == 0 && bet0and2 > 0 && bet1and2 > 0 {
                    ((to0, to1), to2)
                } else if bet0and1 > 0 && bet0and2 == 0 && bet1and2 > 0 {
                    ((to0, to2), to1)
                } else if bet0and1 > 0 && bet0and2 > 0 && bet1and2 == 0 {
                    ((to1, to2), to0)
                } else {
                    return None;
                };
                Some(Bubble {
                    branches,
                    shoot,
                    root,
                })
            })
            .collect()
    }
    fn remove_nodes(&mut self, to_be_removed: &[bool]) {
        let mut next_index = vec![];
        {
            let mut index = 0;
            for &b in to_be_removed.iter() {
                next_index.push(index);
                index += !b as usize;
            }
        }
        // eprint!("{}->", self.indexer.len());
        self.indexer.retain(|_, x| {
            let old_idx = *x;
            *x = next_index[*x];
            !to_be_removed[old_idx]
        });
        // eprintln!("{}", self.indexer.len());
        // eprint!("{}->", self.nodes.len());
        let mut index = 0;
        self.nodes.retain(|_| {
            index += 1;
            !to_be_removed[index - 1]
        });
        // eprintln!("{}", self.nodes.len());
        self.nodes.iter_mut().for_each(|n| {
            n.edges.iter_mut().for_each(|x| x.to = next_index[x.to]);
        });
    }
    fn search_nearest_bubble(
        &self,
        root: usize,
        shoot: usize,
        bubbles: &[Option<Bubble>],
    ) -> Result<usize, ReachedNodeType> {
        let (mut prev, mut current) = (shoot, root);
        while bubbles[current].is_none() {
            let next_nodes = &self.nodes[current].edges;
            if next_nodes.len() == 1 {
                // Fail to find pair. But this is bubble is terminating.
                return Err(ReachedNodeType::Root);
            } else if next_nodes.len() == 2 {
                // Proceed.
                let next = next_nodes.iter().find(|e| e.to != prev).unwrap().to;
                prev = current;
                current = next;
            } else {
                // Fail to find pair. We've reached very complex bubble.
                return Err(ReachedNodeType::ComplexNode);
            }
        }
        if bubbles[current].unwrap().root == prev {
            // We entered current node from the root node. Well done!
            Ok(current)
        } else {
            // We entered current node from either of branches.
            // let (b0, b1) = bubbles[current].unwrap().branches;
            // assert!(prev == b0 || prev == b1);
            Err(ReachedNodeType::BubbleBranch(current))
        }
    }
    /// Resolve crossings in this graph.
    /// The topology of this graph would be changed, and some
    /// of the nodes would be removed.
    pub fn resolve_crossings<T: IntoDeBruijnNodes>(&mut self, reads: &[T]) {
        let mut queue: std::collections::VecDeque<_> = (0..self.nodes.len()).collect();
        let mut to_be_removed = vec![false; self.nodes.len()];
        while let Some(idx) = queue.pop_back() {
            if self.nodes[idx].edges.len() < 4 {
                continue;
            }
            let merged_nodes = match self.resolve_crossing(idx, reads) {
                Some(res) => res,
                None => continue,
            };
            // Success resolving.
            {
                let children: Vec<_> = self.nodes[idx].edges.iter().map(|e| e.to).collect();
                for child in children {
                    self.nodes[child].remove_edge(idx);
                }
            }
            self.nodes[idx].edges.clear();
            to_be_removed[idx] = true;
            for (contract, members) in merged_nodes {
                for member in members {
                    self.nodes[member].push(contract);
                    self.nodes[contract].push(member);
                }
                queue.push_front(contract);
            }
        }
        // eprintln!("{:?}", to_be_removed);
        self.remove_nodes(&to_be_removed);
    }
    fn resolve_crossing<T: IntoDeBruijnNodes>(
        &self,
        center: usize,
        reads: &[T],
    ) -> Option<Vec<(usize, Vec<usize>)>> {
        let neighbors: HashMap<_, _> = self.nodes[center]
            .edges
            .iter()
            .map(|e| e.to)
            .enumerate()
            .map(|(idx, orig)| (orig, idx))
            .collect();
        let mut counts: HashMap<(usize, usize), u32> = HashMap::new();
        for read in reads.iter() {
            let mut occs = vec![false; neighbors.len()];
            for node in read.into_de_bruijn_nodes(self.k) {
                if let Some(&res) = self.indexer.get(&node).and_then(|idx| neighbors.get(idx)) {
                    occs[res] = true;
                }
            }
            for i in 0..neighbors.len() {
                for j in (i + 1)..neighbors.len() {
                    if occs[i] && occs[j] {
                        *counts.entry((i, j)).or_default() += 1;
                    }
                }
            }
        }
        let edges: Vec<_> = counts.into_iter().map(|((x, y), z)| (x, y, z)).collect();
        let (_, used_edges) = maximum_weight_matching(&edges, neighbors.len());
        if used_edges.len() <= 1 {
            return None;
        }
        let mut used_node = vec![false; neighbors.len()];
        for &(n0, n1) in used_edges.iter() {
            used_node[n0] = true;
            used_node[n1] = true;
        }
        let mut cluseters: Vec<_> = used_edges.into_iter().map(|x| vec![x]).collect();
        for (node, _) in used_node.iter().enumerate().filter(|x| !x.1) {
            let (to, _) = match edges
                .iter()
                .filter_map(|&(x, y, w)| {
                    if x == node {
                        Some((y, w))
                    } else if y == node {
                        Some((x, w))
                    } else {
                        None
                    }
                })
                .max_by_key(|x| x.1)
            {
                Some(res) => res,
                None => continue,
            };
            cluseters
                .iter_mut()
                .find(|cl| cl.iter().any(|&(x, y)| x == to || y == to))
                .unwrap()
                .push((to, node));
        }
        // Reverse index.
        let to_original_node_index: HashMap<_, _> = neighbors
            .into_iter()
            .map(|(orig, idx)| (idx, orig))
            .collect();
        Some(
            cluseters
                .into_iter()
                .map(|cl| {
                    let mut degrees: HashMap<usize, u32> = HashMap::new();
                    for (n0, n1) in cl {
                        *degrees.entry(n0).or_default() += 1;
                        *degrees.entry(n1).or_default() += 1;
                    }
                    let (&argmax, _) = degrees.iter().max_by_key(|x| x.1).unwrap();
                    let contract = to_original_node_index[&argmax];
                    degrees.remove(&argmax).unwrap();
                    let members: Vec<_> = degrees
                        .keys()
                        .map(|idx| to_original_node_index[idx])
                        .collect();
                    (contract, members)
                })
                .collect(),
        )
    }
}

#[allow(dead_code)]
// return the median and the median of the absolute error.
fn median_and_mad(xs: &[u64]) -> Option<(u64, u64)> {
    let median = select_nth_by(xs, xs.len() / 2, |&x| x)?;
    let mad = select_nth_by(xs, xs.len() / 2, |&x| {
        if x > median {
            x - median
        } else {
            median - x
        }
    })?;
    Some((median, mad))
}

use std::cmp::{PartialEq, PartialOrd};
#[allow(dead_code)]
fn select_nth_by<T: Clone, F: Fn(&T) -> K, K>(xs: &[T], n: usize, f: F) -> Option<K>
where
    K: PartialOrd + PartialEq,
{
    if xs.len() <= n {
        return None;
    }
    let pivot = f(&xs[xs.len() / 2]);
    let small = xs.iter().filter(|x| f(x) < pivot).count();
    let same = xs.iter().filter(|x| f(x) == pivot).count();
    // Recursive call.
    if n < small {
        // We can remove elements more than `pivot` from `xs`.
        let xs: Vec<_> = xs.iter().filter(|x| f(&x) < pivot).cloned().collect();
        select_nth_by(&xs, n, f)
    } else if small + same <= n {
        let xs: Vec<_> = xs.iter().filter(|x| f(&x) > pivot).cloned().collect();
        select_nth_by(&xs, n - small - same, f)
    } else {
        assert!(small <= n && n < small + same);
        Some(pivot)
    }
}

fn maximum_weight_matching(
    edges: &[(usize, usize, u32)],
    nodes: usize,
) -> (u32, Vec<(usize, usize)>) {
    let mut status = vec![0; edges.len()];
    let (mut max, mut argmax) = (0, vec![]);
    let mut available = vec![true; nodes];
    let mut score = 0;
    let mut stack = vec![0];
    while !stack.is_empty() {
        let current = *stack.last().unwrap();
        if current == edges.len() {
            if max <= score {
                max = score;
                argmax = status.clone();
            }
        } else if status[current] == 0 {
            let (n0, n1, w) = edges[current];
            if available[n0] && available[n1] {
                status[current] = 1;
                stack.push(current + 1);
                score += w;
                available[n0] = false;
                available[n1] = false;
                continue;
            } else {
                status[current] = 2;
                stack.push(current + 1);
                continue;
            }
        } else if status[current] == 1 {
            let (n0, n1, w) = edges[current];
            assert!(!available[n0] && !available[n1]);
            score -= w;
            available[n0] = true;
            available[n1] = true;
            status[current] = 2;
            stack.push(current + 1);
            continue;
        } else {
            assert_eq!(status[current], 2);
            status[current] = 0;
        }
        stack.pop();
    }
    assert_eq!(argmax.len(), edges.len());
    let edges: Vec<_> = edges
        .iter()
        .zip(argmax)
        .filter_map(|(&(n0, n1, _), status)| if status == 1 { Some((n0, n1)) } else { None })
        .collect();
    (max, edges)
}

#[derive(Debug, Clone, Copy)]
struct Bubble {
    // Index of bubble.
    branches: (usize, usize),
    shoot: usize,
    // Root. Where this bubble collapse.
    root: usize,
}

#[derive(Debug, Clone, Copy)]
enum ReachedNodeType {
    Root,
    BubbleBranch(usize),
    ComplexNode,
}

#[cfg(test)]
mod tests {
    use super::*;
    impl AsDeBruijnNode for (u64, u64) {
        fn as_node(w: &[(u64, u64)]) -> Node {
            let first = w.first().unwrap();
            let last = w.last().unwrap();
            let kmer: Vec<_> = if first < last {
                w.iter().copied().collect()
            } else {
                w.iter().rev().copied().collect()
            };
            Node::new(kmer)
        }
    }
    impl IntoDeBruijnNodes for Vec<(u64, u64)> {
        fn into_de_bruijn_nodes(&self, k: usize) -> Vec<Node> {
            self.windows(k).map(AsDeBruijnNode::as_node).collect()
        }
    }
    #[derive(Clone, Copy, Debug)]
    struct TestConfig {
        cl: usize,
        num: usize,
        fail: f64,
        skip: f64,
        max_len: usize,
        min_len: usize,
        unit_len: usize,
    }
    use rand::Rng;
    #[allow(dead_code)]
    fn gen_dataset<R: Rng>(r: &mut R, conf: TestConfig) -> (Vec<Vec<(u64, u64)>>, Vec<usize>) {
        let TestConfig {
            cl,
            num,
            fail,
            skip,
            max_len,
            min_len,
            unit_len,
        } = conf;
        let mut answer = vec![];
        let mut reads = vec![];
        for i in 0..num {
            let cluster = (i % cl) as u64;
            let len = r.gen::<usize>() % (max_len - min_len) + min_len;
            let start = r.gen::<usize>() % (unit_len - len);
            let units: Vec<_> = if r.gen_bool(0.5) {
                (start..=start + len).collect()
            } else {
                let start = start + len;
                (start - len..=start).rev().collect()
            };
            let mut read = vec![];
            for unit in units {
                if r.gen_bool(skip) {
                    continue;
                } else if r.gen_bool(fail) {
                    let cluster = r.gen::<u64>() % cl as u64;
                    read.push((unit as u64, cluster));
                } else {
                    read.push((unit as u64, cluster));
                }
            }
            answer.push(cluster as usize);
            reads.push(read);
        }
        (reads, answer)
    }
    #[test]
    fn construction_test() {
        let read: Vec<Vec<(u64, u64)>> = vec![vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]];
        let graph = DeBruijnGraph::from(&read, 3);
        assert_eq!(graph.nodes.len(), 3, "{:?}", graph);
        let read: Vec<Vec<(u64, u64)>> = vec![
            vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)],
            vec![(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)],
        ];
        let graph = DeBruijnGraph::from(&read, 3);
        assert_eq!(graph.nodes.len(), 3, "{:?}", graph);
        let read: Vec<Vec<(u64, u64)>> = vec![
            vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)],
            vec![(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)],
            vec![(0, 3), (0, 4), (0, 5), (0, 6), (0, 7)],
        ];
        let graph = DeBruijnGraph::from(&read, 3);
        assert_eq!(graph.nodes.len(), 5, "{:?}", graph);
        let read: Vec<Vec<(u64, u64)>> = vec![
            vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)],
            vec![(0, 1), (0, 2), (0, 3), (1, 4), (0, 5), (0, 6), (0, 7)],
            vec![(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)],
            vec![(0, 3), (0, 4), (0, 5), (0, 6), (0, 7)],
        ];
        let graph = DeBruijnGraph::from(&read, 3);
        assert_eq!(graph.nodes.len(), 8, "{:?}", graph);
    }
    #[test]
    fn clustering_test() {
        let read: Vec<Vec<(u64, u64)>> = vec![vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]];
        let graph = DeBruijnGraph::from(&read, 3);
        assert_eq!(graph.clustering(0).len(), 1, "{:?}", graph);
        let read: Vec<Vec<(u64, u64)>> = vec![
            vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)],
            vec![(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)],
        ];
        let graph = DeBruijnGraph::from(&read, 3);
        assert_eq!(graph.clustering(0).len(), 1, "{:?}", graph);
        let read: Vec<Vec<(u64, u64)>> = vec![
            vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)],
            vec![(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)],
            vec![(0, 3), (0, 4), (0, 5), (0, 6), (0, 7)],
        ];
        let graph = DeBruijnGraph::from(&read, 3);
        assert_eq!(graph.clustering(0).len(), 1, "{:?}", graph);
        let read: Vec<Vec<(u64, u64)>> = vec![
            vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)],
            vec![(0, 1), (0, 2), (0, 3), (1, 4), (0, 5), (0, 6), (0, 7)],
            vec![(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)],
            vec![(0, 3), (0, 4), (0, 5), (0, 6), (0, 7)],
        ];
        let graph = DeBruijnGraph::from(&read, 3);
        assert_eq!(graph.clustering(0).len(), 1, "{:?}", graph);
    }
    #[test]
    fn clustering_test_2() {
        // Case1
        let reads: Vec<Vec<(u64, u64)>> = vec![
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 0), (3, 1), (4, 1), (5, 1)],
            vec![(4, 1), (3, 1), (2, 1), (1, 1)],
        ];
        let graph = DeBruijnGraph::from(&reads, 3);
        assert_eq!(graph.clustering(0).len(), 2, "\n{:?}", graph);
        // Case 2
        let reads: Vec<Vec<(u64, u64)>> = vec![
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0)],
            vec![(0, 1), (1, 1), (2, 0), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 0), (3, 1), (4, 1), (5, 1)],
            vec![(4, 1), (3, 1), (2, 0), (1, 1)],
        ];
        let graph = DeBruijnGraph::from(&reads, 3);
        assert_eq!(graph.clustering(0).len(), 2, "\n{:?}", graph);
        // Case 3
        let reads: Vec<Vec<(u64, u64)>> = vec![
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(5, 1), (4, 1), (3, 1), (2, 1), (1, 1), (0, 1)],
            vec![(0, 1), (1, 1), (2, 2), (3, 1), (4, 1), (5, 1)],
            vec![(4, 1), (3, 1), (2, 2), (1, 1)],
        ];
        let graph = DeBruijnGraph::from(&reads, 3);
        assert_eq!(graph.clustering(0).len(), 2, "\n{:?}", graph);
        // Case 4
        let reads: Vec<Vec<(u64, u64)>> = vec![
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0), (1, 0), (0, 0)],
            vec![(0, 0), (1, 0), (3, 0), (4, 0), (5, 0)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
        ];
        let graph = DeBruijnGraph::from(&reads, 3);
        assert_eq!(graph.clustering(0).len(), 2, "\n{:?}", graph);
    }
    #[test]
    fn clustering_test_3() {
        let reads: Vec<Vec<(u64, u64)>> = vec![
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0), (1, 0), (0, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0), (1, 0), (0, 0)],
            vec![(1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(5, 1), (4, 1), (3, 1), (2, 1), (1, 1), (1, 0)],
            vec![(5, 1), (4, 1), (3, 1), (2, 1), (1, 1), (1, 0)],
            vec![(4, 1), (3, 1), (2, 1), (1, 1), (1, 0)],
            vec![(1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 0), (3, 1), (4, 1), (5, 1)],
            vec![(5, 1), (4, 1), (3, 1), (2, 0), (1, 1), (0, 1)],
        ];
        let graph = DeBruijnGraph::from(&reads, 3).clean_up(1);
        assert_eq!(graph.clustering(0).len(), 2, "\n{:?}", graph);
        let reads: Vec<Vec<(u64, u64)>> = vec![
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0), (1, 0), (0, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0), (1, 0), (0, 0)],
            vec![(1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(5, 1), (4, 1), (3, 1), (2, 1), (1, 1), (0, 1)],
            vec![(5, 1), (4, 1), (3, 1), (2, 1), (1, 1), (0, 1)],
            vec![(4, 1), (3, 1), (2, 1), (1, 1), (0, 1)],
            vec![(1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 0), (3, 0), (4, 0), (5, 1)],
            vec![(5, 1), (4, 0), (3, 0), (2, 0), (1, 1), (0, 1)],
        ];
        let graph = DeBruijnGraph::from(&reads, 3);
        let graph = graph.clean_up(1);
        assert_eq!(graph.clustering(0).len(), 2, "\n{:?}", graph);
    }
    #[test]
    fn path_clustering_test_large_noisy() {
        use rand::SeedableRng;
        use rand_xoshiro::Xoshiro256StarStar;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3205);
        let conf = TestConfig {
            cl: 2,
            num: 2000,
            fail: 0.02,
            skip: 0.02,
            max_len: 40,
            min_len: 5,
            unit_len: 500,
        };
        let (reads, answer) = gen_dataset(&mut rng, conf);
        let graph = DeBruijnGraph::from(&reads, 3);
        let graph = graph.clean_up_auto();
        eprintln!("{:?}", graph);
        let mut components = graph.clustering(1);
        components.retain(|c| c.len() > 100);
        assert_eq!(components.len(), 2);
        let preds: Vec<_> = reads
            .iter()
            .filter_map(|read| {
                components
                    .iter()
                    .map(|c| read.iter().filter(|x| c.contains(x)).count())
                    .enumerate()
                    .max_by_key(|x| x.1)
            })
            .map(|x| x.0)
            .collect();
        assert_eq!(preds.len(), answer.len());
        let correct = answer
            .iter()
            .zip(preds.iter())
            .filter(|&(ans, pred)| ans == pred)
            .count();
        eprintln!("{}/{}", correct, reads.len());
        for ((read, assign), ans) in reads.iter().zip(preds.iter()).zip(answer.iter()) {
            eprintln!("{}\t{}\t{}", ans, assign, read.len());
        }
        let correct = correct.max(reads.len() - correct);
        assert!(correct > reads.len() * 8 / 10);
    }
    #[test]
    fn path_clustering_test_large_hard() {
        use rand::SeedableRng;
        use rand_xoshiro::Xoshiro256StarStar;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3205);
        let conf = TestConfig {
            cl: 2,
            num: 5000,
            fail: 0.02,
            skip: 0.02,
            max_len: 40,
            min_len: 10,
            unit_len: 800,
        };
        let (reads, answer) = gen_dataset(&mut rng, conf);
        let graph = DeBruijnGraph::from(&reads, 3);
        let graph = graph.clean_up_auto();
        eprintln!("{:?}", graph);
        let mut components = graph.clustering(1);
        components.retain(|cl| cl.len() > 100);
        assert_eq!(components.len(), 2);
        let preds: Vec<_> = reads
            .iter()
            .filter_map(|read| {
                components
                    .iter()
                    .map(|c| read.iter().filter(|x| c.contains(x)).count())
                    .enumerate()
                    .max_by_key(|x| x.1)
            })
            .map(|x| x.0)
            .collect();
        assert_eq!(preds.len(), answer.len());
        let correct = answer
            .iter()
            .zip(preds.iter())
            .filter(|&(ans, pred)| ans == pred)
            .count();
        eprintln!("{}/{}", correct, reads.len());
        for ((read, assign), ans) in reads.iter().zip(preds.iter()).zip(answer.iter()) {
            eprintln!("{}\t{}\t{}", ans, assign, read.len());
        }
        let correct = correct.max(reads.len() - correct);
        assert!(correct > reads.len() * 8 / 10);
    }
    #[test]
    fn path_clustering_test_short() {
        use rand::SeedableRng;
        use rand_xoshiro::Xoshiro256StarStar;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3205);
        let conf = TestConfig {
            cl: 2,
            num: 200,
            fail: 0.01,
            skip: 0.01,
            max_len: 20,
            min_len: 10,
            unit_len: 50,
        };
        let (reads, answer) = gen_dataset(&mut rng, conf);
        let graph = DeBruijnGraph::from(&reads, 3);
        let graph = graph.clean_up_auto();
        eprintln!("{:?}", graph);
        let mut components = graph.clustering(0);
        components.retain(|c| c.len() > 20);
        assert_eq!(components.len(), 2);
        let preds: Vec<_> = reads
            .iter()
            .filter_map(|read| {
                components
                    .iter()
                    .map(|c| read.iter().filter(|x| c.contains(x)).count())
                    .enumerate()
                    .max_by_key(|x| x.1)
            })
            .map(|x| x.0)
            .collect();
        assert_eq!(preds.len(), answer.len());
        let correct = answer
            .iter()
            .zip(preds.iter())
            .filter(|&(ans, pred)| ans == pred)
            .count();
        eprintln!("{}/{}", correct, reads.len());
        for ((read, assign), ans) in reads.iter().zip(preds.iter()).zip(answer.iter()) {
            eprintln!("{}\t{}\t{}", ans, assign, read.len());
        }
        let correct = correct.max(reads.len() - correct);
        assert!(correct > reads.len() * 8 / 10);
    }
    // Bubble resolving functionality.
    #[test]
    fn bubble_works() {
        let reads = vec![vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]];
        let mut graph = DeBruijnGraph::from(&reads, 3);
        graph.resolve_bubbles(&reads)
    }
    fn sample_from_to(
        template: &[(u64, u64)],
        from: usize,
        to: usize,
        rev: bool,
    ) -> Vec<(u64, u64)> {
        if rev {
            template[from..to].iter().rev().copied().collect()
        } else {
            template[from..to].to_vec()
        }
    }
    #[test]
    fn bubble_case_0() {
        let templates: Vec<_> = vec![vec![0; 8], vec![1, 1, 0, 0, 0, 0, 1, 0]];
        let templates: Vec<Vec<(u64, u64)>> = templates
            .iter()
            .map(|ts| {
                ts.iter()
                    .enumerate()
                    .map(|(idx, &x)| (idx as u64, x as u64))
                    .collect()
            })
            .collect();
        let reads = vec![
            sample_from_to(&templates[0], 0, 8, false),
            sample_from_to(&templates[1], 0, 8, false),
            sample_from_to(&templates[0], 0, 8, true),
            sample_from_to(&templates[1], 0, 8, true),
        ];
        let mut graph = DeBruijnGraph::from(&reads, 3);
        eprintln!("{:?}", graph);
        let bubbles = graph.enumerate_bubbles(&reads);
        let num_bubbles = bubbles.iter().filter(|e| e.is_some()).count();
        assert_eq!(num_bubbles, 2, "{:?}", graph);
        let pair_idx = graph.search_nearest_bubble(3, 2, &bubbles).unwrap();
        let bubble = bubbles[pair_idx].unwrap();
        assert_eq!(bubble.shoot, 3);
        graph.resolve_bubbles(&reads);
        assert_eq!(graph.clustering(0).len(), 2, "{:?}", graph);
        graph.validate();
    }
    #[test]
    fn bubble_case_1() {
        let templates: Vec<_> = vec![
            vec![0; 9],
            vec![0, 1, 0, 0, 0, 0, 1, 1, 0],
            vec![0, 1, 0, 0, 0, 0, 1, 1, 2],
        ];
        let templates: Vec<Vec<(u64, u64)>> = templates
            .iter()
            .map(|ts| {
                ts.iter()
                    .enumerate()
                    .map(|(idx, &x)| (idx as u64, x as u64))
                    .collect()
            })
            .collect();
        let reads = vec![
            sample_from_to(&templates[0], 0, 9, false),
            sample_from_to(&templates[1], 0, 9, false),
            sample_from_to(&templates[2], 0, 9, false),
            sample_from_to(&templates[0], 0, 9, true),
            sample_from_to(&templates[1], 0, 9, true),
            sample_from_to(&templates[2], 0, 9, true),
        ];
        let mut graph = DeBruijnGraph::from(&reads, 3);
        eprintln!("{:?}", graph);
        let bubbles = graph.enumerate_bubbles(&reads);
        let num_bubbles = bubbles.iter().filter(|e| e.is_some()).count();
        assert_eq!(num_bubbles, 3, "{:?}", graph);
        graph.resolve_bubbles(&reads);
        assert_eq!(graph.clustering(0).len(), 3, "{:?}", graph);
        graph.validate();
    }
    #[test]
    fn bubble_case_2() {
        let template_a: Vec<(u64, u64)> = (0..11).map(|idx| (idx, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let template_c: Vec<(u64, u64)> = vec![1, 1, 0, 0, 1, 1, 2, 1]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64 + 3, c))
            .collect();
        let reads = vec![
            sample_from_to(&template_a, 0, 11, false),
            sample_from_to(&template_a, 0, 11, true),
            sample_from_to(&template_b, 0, 11, false),
            sample_from_to(&template_b, 0, 11, true),
            sample_from_to(&template_c, 0, 8, false),
            sample_from_to(&template_c, 0, 8, true),
        ];
        let mut graph = DeBruijnGraph::from(&reads, 3);
        eprintln!("{:?}", graph);
        let bubbles = graph.enumerate_bubbles(&reads);
        let num_bubbles = bubbles.iter().filter(|e| e.is_some()).count();
        assert_eq!(num_bubbles, 4, "{:?}", graph);
        graph.resolve_bubbles(&reads);
        assert_eq!(graph.clustering(0).len(), 3, "{:?}", graph);
        graph.validate();
    }
    #[test]
    fn bubble_case_3() {
        let template_a: Vec<(u64, u64)> = (0..13).map(|idx| (idx, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let reads = vec![
            sample_from_to(&template_a, 0, 13, false),
            sample_from_to(&template_a, 0, 13, true),
            sample_from_to(&template_b, 0, 13, false),
            sample_from_to(&template_b, 0, 13, true),
        ];
        let mut graph = DeBruijnGraph::from(&reads, 3);
        eprintln!("{:?}", graph);
        let bubbles = graph.enumerate_bubbles(&reads);
        let num_bubbles = bubbles.iter().filter(|e| e.is_some()).count();
        assert_eq!(num_bubbles, 4, "{:?}", graph);
        graph.resolve_bubbles(&reads);
        assert_eq!(graph.clustering(0).len(), 2, "{:?}", graph);
        graph.validate();
    }
    #[test]
    fn bubble_case_4() {
        let template_a: Vec<(u64, u64)> = (0..12).map(|idx| (idx, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let reads = vec![
            sample_from_to(&template_a, 0, 12, false),
            sample_from_to(&template_a, 0, 12, true),
            sample_from_to(&template_b, 0, 12, false),
            sample_from_to(&template_b, 0, 12, true),
        ];
        let mut graph = DeBruijnGraph::from(&reads, 3);
        eprintln!("{:?}", graph);
        let bubbles = graph.enumerate_bubbles(&reads);
        let num_bubbles = bubbles.iter().filter(|e| e.is_some()).count();
        assert_eq!(num_bubbles, 3, "{:?}", graph);
        graph.resolve_bubbles(&reads);
        assert_eq!(graph.clustering(0).len(), 2, "{:?}", graph);
        graph.validate();
    }
    #[test]
    fn bubble_case_5() {
        let template_a: Vec<(u64, u64)> = (0..9).map(|idx| (idx + 2, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let template_c: Vec<(u64, u64)> = vec![2, 2, 1, 1, 0, 0, 0, 0, 0, 1, 1, 2, 1]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let reads = vec![
            sample_from_to(&template_a, 0, 9, false),
            sample_from_to(&template_a, 0, 9, true),
            sample_from_to(&template_b, 0, 13, false),
            sample_from_to(&template_b, 0, 13, true),
            sample_from_to(&template_c, 0, 13, false),
            sample_from_to(&template_c, 0, 13, true),
        ];
        let mut graph = DeBruijnGraph::from(&reads, 3);
        eprintln!("{:?}", graph);
        let bubbles = graph.enumerate_bubbles(&reads);
        let num_bubbles = bubbles.iter().filter(|e| e.is_some()).count();
        assert_eq!(num_bubbles, 4, "{:?}", graph);
        graph.resolve_bubbles(&reads);
        assert_eq!(graph.clustering(0).len(), 3, "{:?}", graph);
        graph.validate();
    }
    #[test]
    fn maximum_weight_matching_test() {
        let nodes = 4;
        let edges = vec![(0, 1, 2), (2, 3, 2)];
        let (score, result) = maximum_weight_matching(&edges, nodes);
        assert_eq!(result.len(), 2);
        assert_eq!(score, 4);
        let nodes = 4;
        let edges = vec![
            (0, 1, 1),
            (1, 2, 1),
            (2, 3, 1),
            (3, 0, 1),
            (0, 2, 3),
            (1, 3, 3),
        ];
        let (score, result) = maximum_weight_matching(&edges, nodes);
        assert_eq!(result.len(), 2);
        assert_eq!(score, 6);
        let nodes = 6;
        let edges = vec![
            (0, 1, 4),
            (0, 2, 2),
            (1, 2, 2),
            (2, 3, 8),
            (3, 4, 2),
            (4, 5, 5),
            (3, 5, 2),
        ];
        let (score, result) = maximum_weight_matching(&edges, nodes);
        assert_eq!(result.len(), 3);
        assert_eq!(score, 17);
    }
    #[test]
    fn resolve_crossing_0() {
        let template_a: Vec<(u64, u64)> = (0..7).map(|idx| (idx, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![0, 1, 0, 0, 0, 1, 0]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let reads = vec![
            sample_from_to(&template_a, 0, 7, false),
            sample_from_to(&template_a, 0, 7, true),
            sample_from_to(&template_b, 0, 7, false),
            sample_from_to(&template_b, 0, 7, true),
        ];
        let mut graph = DeBruijnGraph::from(&reads, 3);
        eprintln!("{:?}", graph);
        graph.resolve_crossings(&reads);
        assert_eq!(graph.clustering(0).len(), 2, "{:?}", graph);
        graph.validate();
    }
    #[test]
    fn resolve_crossing_1() {
        let template_a: Vec<(u64, u64)> = (0..7).map(|idx| (idx, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![0, 1, 0, 0, 0, 1, 0]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let template_c: Vec<(u64, u64)> = vec![0, 2, 0, 0, 0, 2, 0]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let reads = vec![
            sample_from_to(&template_a, 0, 7, false),
            sample_from_to(&template_a, 0, 7, true),
            sample_from_to(&template_b, 0, 7, false),
            sample_from_to(&template_b, 0, 7, true),
            sample_from_to(&template_c, 0, 7, false),
            sample_from_to(&template_c, 0, 7, true),
        ];
        let mut graph = DeBruijnGraph::from(&reads, 3);
        eprintln!("{:?}", graph);
        graph.resolve_crossings(&reads);
        assert_eq!(graph.clustering(0).len(), 3, "{:?}", graph);
        graph.validate();
    }
    #[test]
    fn resolve_crossing_2() {
        let template_a: Vec<(u64, u64)> = (0..7).map(|idx| (idx, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![0, 1, 0, 0, 0, 1, 0]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let template_c: Vec<(u64, u64)> = vec![0, 1, 0, 0, 0, 2, 0]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let reads = vec![
            sample_from_to(&template_a, 0, 7, false),
            sample_from_to(&template_a, 0, 7, true),
            sample_from_to(&template_b, 0, 7, false),
            sample_from_to(&template_b, 0, 7, true),
            sample_from_to(&template_c, 0, 7, false),
            sample_from_to(&template_c, 0, 7, true),
        ];
        let mut graph = DeBruijnGraph::from(&reads, 3);
        eprintln!("{:?}", graph);
        graph.resolve_crossings(&reads);
        assert_eq!(graph.clustering(0).len(), 2, "{:?}", graph);
        graph.validate();
    }
    #[test]
    fn resolve_crossing_3() {
        let template_a: Vec<(u64, u64)> = (0..7).map(|idx| (idx, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![0, 1, 0, 0, 0, 1, 0]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let template_c: Vec<(u64, u64)> = vec![0, 2, 0, 0, 0, 2, 0]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let templates: Vec<Vec<(u64, u64)>> = vec![
            vec![0, 0, 0, 0, 0, 1, 0],
            vec![0, 2, 0, 0, 0, 1, 0],
            vec![0, 0, 0, 0, 0, 2, 0],
            vec![0, 1, 0, 0, 0, 0, 0],
        ]
        .into_iter()
        .map(|template| {
            template
                .into_iter()
                .enumerate()
                .map(|(idx, c)| (idx as u64, c))
                .collect()
        })
        .collect();
        let mut reads = vec![
            sample_from_to(&template_a, 0, 7, false),
            sample_from_to(&template_a, 0, 7, true),
            sample_from_to(&template_b, 0, 7, false),
            sample_from_to(&template_b, 0, 7, true),
            sample_from_to(&template_c, 0, 7, false),
            sample_from_to(&template_c, 0, 7, true),
            sample_from_to(&template_a, 0, 7, false),
            sample_from_to(&template_a, 0, 7, true),
            sample_from_to(&template_b, 0, 7, false),
            sample_from_to(&template_b, 0, 7, true),
            sample_from_to(&template_c, 0, 7, false),
            sample_from_to(&template_c, 0, 7, true),
        ];
        for template in templates {
            reads.push(sample_from_to(&template, 0, 7, true));
            reads.push(sample_from_to(&template, 0, 7, false));
        }
        let mut graph = DeBruijnGraph::from(&reads, 3);
        eprintln!("{:?}", graph);
        graph.resolve_crossings(&reads);
        assert_eq!(graph.clustering(0).len(), 3, "{:?}", graph);
        graph.validate();
    }
    #[test]
    fn resolve_crossing_4() {
        let template_a: Vec<(u64, u64)> = (0..7).map(|idx| (idx, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![0, 1, 0, 0, 0, 0, 1]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let template_c: Vec<(u64, u64)> = vec![0, 2, 0, 0, 0, 2, 0]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let reads = vec![
            sample_from_to(&template_a, 0, 7, false),
            sample_from_to(&template_a, 0, 7, true),
            sample_from_to(&template_b, 0, 7, false),
            sample_from_to(&template_b, 0, 7, true),
            sample_from_to(&template_c, 0, 7, false),
            sample_from_to(&template_c, 0, 7, true),
        ];
        let mut graph = DeBruijnGraph::from(&reads, 3);
        eprintln!("{:?}", graph);
        graph.resolve_crossings(&reads);
        assert_eq!(graph.clustering(0).len(), 3, "{:?}", graph);
        graph.validate();
    }
    #[test]
    fn bubble_resolve_and_crossing_resolve() {
        let template_a: Vec<(u64, u64)> = (0..12).map(|idx| (idx, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let reads = vec![
            sample_from_to(&template_a, 0, 11, false),
            sample_from_to(&template_a, 0, 11, true),
            sample_from_to(&template_b, 0, 11, false),
            sample_from_to(&template_b, 0, 11, true),
        ];
        let mut graph = DeBruijnGraph::from(&reads, 3);
        eprintln!("Before:{:?}", graph);
        graph.resolve_crossings(&reads);
        graph.resolve_bubbles(&reads);
        assert_eq!(graph.clustering(0).len(), 2, "{:?}", graph);
        graph.validate();
    }
}
