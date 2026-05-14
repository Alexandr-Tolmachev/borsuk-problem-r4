# Reducing the upper bound for the Borsuk number in $\mathbb{R}^4$ to 8
Alexander Tolmachev, Vsevolod Voronov

<p align="center">
  <a href="https://arxiv.org/abs/xxxx.xxxxx"><b>📄 Paper</b></a>
</p>

---
> 🧠 **Abstract:**
> The Borsuk number $b(n)$ of $n$-dimensional Euclidean space $\mathbb{R}^n$ is the smallest integer such that any set $F \subset \mathbb{R}^n$ of unit diameter can be partitioned into $b(n)$ subsets of strictly smaller diameter. For $n=4$, the best known upper bound $b(4) \leq 9$ follows from a construction by M. Lassak (1982). In the present paper, we construct partitions of several variants of the truncated Lassak cover into 8 parts of diameter less than 1, thereby showing that $b(4) \leq 8$. 
---

#### 🔽 Installation

First, clone the repository and install the required Python libraries. Optionally (but recommended), create a dedicated virtual environment to avoid dependency conflicts:

```bash
git clone https://github.com/Alexandr-Tolmachev/borsuk-problem-r4.git
cd borsuk-problem-r4

# Optional but recommended
python -m venv .venv
source .venv/bin/activate  # Linux/macOS
# .venv\Scripts\activate   # Windows

pip install -r requirements.txt
```

#### 🔽 Validation of partitions

The folder `covers` contains four files with detailed descriptions of the partition of each element of the constructed universal covering system (UCS) into 8 parts.

To run the validation procedure, use the following command:

```bash
python validation.py <path_to_file> -g <grid_size>
```

Here, `path_to_file` indicates the path to the file containing the cover description, and `grid_size` is the integer discretization parameter $(g \ge 16)$ that controls the accuracy of the outer polyhedral approximation used during the diameter estimation procedure.

**Example:**

Running the command `validation.py covers/cover1 -g 17` starts the validation procedure for the first cover set with a grid size of 17 for the polyhedral approximation construction.


#### 🔽 Final verification
To validate and reproduce the theoretical result from our [paper](https://arxiv.org/abs/xxxx.xxxxx), you should run this procedure for all covers in the corresponding folder (4 items) with $g=17$:

```bash
python validation.py covers/cover1 -g 17
python validation.py covers/cover2 -g 17
python validation.py covers/cover3 -g 17
python validation.py covers/cover4 -g 17
```

#### 🔽 Visualization

To visualize the 3D projection of the eight parts of the cover partition, we provide code that generates GIF animations. These animations illustrate the projection onto a moving 4D plane, parameterized by its normal vector while varying the offset (shift).

```bash
python validation.py <path_to_file> -g <grid_size> \
                                    -normal <normal_vector> \
                                    -gif_duration <gif_duration> \
                                    -gif_path <output_path> \
                                    -num_frames <num_frames>
```
Here, 
- `path_to_file` indicates the path to the file containing the cover description,

- `normal` defines the coordinates of the normal vector (using this format without additional blanks: [1,2,3,4] or 1,2,3,4), 

- `grid_size` is the integer discretization parameter (use small values $3 \le g \le 6$ for better visualization), 

- `gif_duration` is the video duration (in seconds), 

- `output_path` indicates the path to save final animation, 

- `num_frames` corresponds to the number of frames in the video (number of plane shifts along the direction parameterized via `normal` vector).

For instance, we provide several examples for GIF creation:

```bash
# Basic usage with defaults
python visualization.py covers/cover1

# With custom grid size and GIF parameters
python visualization.py covers/cover1 -g 4 -gif_duration 3 -gif_path output.gif

# Full custom configuration
python visualization.py covers/cover1 -g 4 -normal [1,1,1,1] -gif_duration 6 \
                                      -num_frames 90 -gif_path animation.gif

# with another option for normal vector input
python visualization.py covers/cover1 -g 4 -normal 1,1,1,1 -gif_duration 6 \
                                      -num_frames 180 -gif_path animation.gif
```
**Example:** projections of cover1 partition elements (8 items) into several moving 4D planes:

<div style="display: flex; gap: 10px; justify-content: center;">
  <img src="gif_examples/cover1_proj1.gif" alt="GIF 1" width="250">
  <img src="gif_examples/cover1_proj2.gif" alt="GIF 2" width="250">
  <img src="gif_examples/cover1_proj3.gif" alt="GIF 3" width="250">
</div>

---

## 📚 Citation

If you find our work or code useful, please cite:

```bibtex
@article{tolmachev2026borsukR4,
  title={Breaking a Barrier in the Borsuk Problem: A New Upper Bound for $R^4$},
  author={Tolmachev, Alexander and Voronov, Vsevolod},
  journal={arXiv preprint arXiv:xxxx.xxxxx},
  year={2026}
}
```
