# How to use CVD (the simple version)

*Written so a 13-year-old can follow it. The full details are in the [README](../README.md).*

## What is this thing?

Your body is built from tiny cells, and each cell follows a recipe book called **DNA**. Each
recipe in the book is a **gene**. Sometimes a gene breaks, and that can cause a disease like
cancer.

**Gene therapy** is an idea for fixing that: deliver a *good copy* of a recipe into the cells.
To get it inside, scientists use a **vector**, which is a delivery truck made from a harmless,
stripped-down virus. The truck carries a small package of DNA instead of germs.

**CVD is a design program.** You tell it what you want the delivery truck to carry, and it draws
up the plan: a blueprint listing every part, in order, with a circular map. It does not make
anything real. It only makes the plan, on your computer.

> **Important:** this is a school-and-research design tool. Building real viruses needs a
> professional lab, trained people and official safety approval. Never try that on your own.

## What you need

- A computer with **Python** (version 3.8 or newer). Get it free from python.org.
- The CVD folder, unzipped.

That's all. You don't have to install anything else.

## Opening the program

Open a terminal (a window where you type commands) inside the CVD folder and type:

```
python -m cvd
```

A window opens. If your computer can't open a window, use the question-and-answer version:

```
python -m cvd --wizard
```

(There is also a one-file version: `python cvd.pyz` does the same thing.)

## Making your first design (window version)

1. **Pick a vector.** Think of it as the type of delivery truck. *Lentivirus* and *AAV* are two
   kinds. *Plasmid* is just the plain DNA circle with no virus parts.
2. **Pick a strategy.** This is the job you want done. For example, *Tumor suppressor
   replacement* means "put back a good copy of a gene that cancer broke."
3. **Pick a gene** from the list, like TP53, or type your own gene name.
4. **Look at the map.** The colored circle is your plan. Each colored piece is one part, and
   the table next to it lists every part with its size.
5. **Read the checks.** The program warns you if something looks wrong, such as a package that
   is too big for the truck.
6. **Save your files.** The program writes the plan as several files (a map picture and
   DNA-format files that biologists' software can open).

Tip: parts shown as "N" are placeholders. They are the right *size* but not the real DNA
letters yet. Scientists swap in the verified real ones before anything is built.

## The ranked options tab (the cool part)

Click the **Vector options** tab. The program builds several complete plans for your gene and
**ranks them from best fit to worst fit**, like a leaderboard.

- Each plan has a small circular map. **Click a map** to open the full plan.
- Each plan has a score out of 100. A truck that is too small for the package gets 0.
- Click a plan in the list and the **Why this ranking** tab explains why it is above or below
  every other plan. Nothing is a secret: the reasons come from the actual scores.
- The **Modifications** tab lists every change the program made to improve that plan, such as
  swapping in a shorter part so everything fits. Each change shows why, and what it costs.
- If you want to see the plans *without* the improvements, untick **Apply suggested
  modifications**.
- Like a plan? Press **Use this option in the designer**.

Remember: the scores are a guide made from rules. They compare the plans against each other.
They do not predict what would really happen in a lab.

## Doing it by typing commands instead

```
python -m cvd --list
```
shows every gene and part you can choose from.

```
python -m cvd --strategy rep --gene TP53 --setting "in vivo" --rank --explain 1
```
ranks the options for the TP53 gene and explains why option 1 is on top. ("In vivo" means the
treatment goes inside the body. "Ex vivo" means cells are treated outside the body first.)

```
python -m cvd --gene TP53
```
makes one design and saves the files.

## The chromosome-counting part (aneuploidy)

Healthy human cells have 23 pairs of chromosomes, which are like 23 pairs of big recipe
volumes. Cancer cells often have *extra or missing pieces* of those volumes. That is called
**aneuploidy**, and it can make cancer cells look different from normal ones.

CVD can read a file of those counts, tell you which chromosome arms are gained or lost, give
a score, and suggest designs that help the immune system notice the odd cells. To try it with
made-up practice data:

```
python -m cvd --example-cnv demo.tsv
python -m cvd --cnv demo.tsv
```

## If something goes wrong

- **"python is not recognized":** Python isn't installed, or isn't on your PATH. Reinstall it
  and tick "Add Python to PATH."
- **No window opens:** your Python may be missing tkinter. Use `python -m cvd --wizard`.
- **A message that starts with "Error":** read it. It usually says exactly what to fix, like
  a typo in a gene name.

## Words to know

| Word | Plain meaning |
|---|---|
| DNA | the recipe book inside every cell |
| Gene | one recipe in the book |
| Vector | the delivery truck for the DNA package |
| Promoter | the "on switch" for a gene |
| Insert | the new gene you are delivering |
| Cassette | the package inside the truck |
| Aneuploidy | extra or missing chromosome pieces |

Have fun, and ask a teacher or scientist to check your ideas before taking them further.
