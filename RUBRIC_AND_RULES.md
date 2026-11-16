# Rubric and rules: Veritas AI International High School Hackathon 2026

A working summary of the organisers' documents. The originals are in `docs/` and are the authority if anything here differs.

---

## Key dates (US Central)

| Event | When |
| --- | --- |
| Hacking started | Fri Nov 13, 11:00 PM |
| **Submissions close** | **Sun Nov 15, 10:59 PM** (11:59 PM Eastern) |
| Shortlist and finalists emailed | Wed Nov 18 |
| Live finals and awards | Sat Nov 21, 10:00 AM (11:00 AM Eastern) |

---

## The four judging criteria

Every track is scored on the same four.

| Criterion | What judges look for | How to show it |
| --- | --- | --- |
| **Technical quality** | The project works, the approach is sound, and it was tested fairly | Compare with a baseline, report results on unseen test data, show where it fails |
| **Impact** | It tackles a real problem for real people | Explain who it helps and what would change for them |
| **Creativity** | A fresh approach, not a repeated tutorial | Apply a known technique to an unusual dataset, question or group of users |
| **Clarity** | A judge can understand it in two minutes | A clear demo, a concise write-up, one main takeaway |

Judges score an honest, modest result higher than an unexplained high one. Results that look too good usually mean test data leaked into training.

---

## What every project must show

1. **Who is it for?** One specific user and the moment they would reach for the tool.
2. **What goes in, and what comes out?** A real input producing a result.
3. **Does it beat a simple baseline?** The baseline's score next to the model's score.
4. **Where does it fail?** At least one example the model gets wrong, with a reason.

All results must come from test data the model did not see during training.

---

## AI + Health track

**Challenge:** build an AI tool that helps one specific person catch a health problem sooner, or manage an ongoing one better.

**Required**

- Name the user and explain how they would use the tool.
- Report missed cases (false negatives) and false alarms (false positives) separately, not only accuracy.
- Include this line in the write-up, word for word: "This project is a prototype and is not medical advice."
- Public datasets only. No real patient records.

**Judges pay extra attention to**

- How clearly the user is identified.
- How honestly missed cases and false alarms are reported.
- Whether the model was checked across groups such as age and sex.

---

## The rules

1. **Teams** of one to three high school students. Everyone listed must have contributed.
2. **Build during the weekend.** Learning tools and browsing datasets beforehand is fine; extending an earlier project is not.
3. **AI must be central.** A regular app built with AI coding tools does not count on its own.
4. **Allowed:** public datasets, pretrained models, AI APIs, open-source libraries, tutorial code, coding assistants and chatbots. Cite all of them. The team must understand every part and be ready to explain it.
5. **Not allowed:** submitting someone else's project; collecting personal information without permission.
6. **One project, one track.** Resubmitting is allowed; judges see the latest version.
7. **Disclaimers** for health and finance projects.
8. **Honest results,** including weak ones.

---

## Own idea or suggested project

- Suggested projects are compared with each other and scored on what the team added.
- A team following the steps exactly scores toward the lower end.
- Adding your own user, your own data or your own question scores higher.
- An own idea writes "Our own idea" in the "What we added" section.

---

## What to submit

| Item | Required? | Details |
| --- | --- | --- |
| Code | Yes | A Colab notebook or GitHub repository anyone with the link can view. Notebook saved with outputs visible. |
| Write-up | Yes | About one page, on the template below. A Google Doc, the GitHub README, or an uploaded file. |
| Demo video | Optional, recommended | Up to 2 minutes. YouTube unlisted, or Google Drive with link sharing. |
| Team details | Yes | Track, and each member's name, email and school. |
| Credits | Yes | Every dataset, pretrained model, tutorial and AI tool, with links. |

A submission is complete once it has working links to both the code and the write-up.

**Before submitting**

- Open every link in a private window to check it works without signing in.
- Check every teammate's name and email.
- Submit a first version a few hours early, then resubmit.

### Form fields seen in the sample submission

Project title; building as a team or alone; team name; team members; track; starting point; one-sentence pitch; code link; write-up link; demo video link; what you added; datasets, models and AI tools; new to AI?

---

## Write-up template

```
# Project name
**Track:** Health / Finance / Climate / Languages
**Team:** names

## The problem
Who has this problem, and why does it matter? (2-3 sentences)

## What we built
What goes in, what comes out, and how someone would use it.

## What we added
If you started from a suggested project: which one, and what you added or changed.
Otherwise, write "Our own idea".

## How it works
Data: which dataset, how many examples, how we cleaned it.
Model: what we used and why.

## Results
Our baseline scored X. Our model scored Y on held-out test data.
Where it fails, and why we think so.

## What we'd do next
One or two concrete next steps with more time.

## Sources
Datasets, pretrained models, tutorials and AI tools used, with links.
```

## Video structure (2 minutes)

| Part | Length | What to show |
| --- | --- | --- |
| The problem | About 15 seconds | Who it is for and why it matters |
| The demo | About 60 seconds | A real input going in and the result coming out |
| How it works | About 30 seconds | Data, model, and comparison with the baseline |
| Next steps | About 15 seconds | One limitation and what comes next |

A screen recording with a voice-over is enough.

---

## How judging runs

- **Round 1:** the Veritas AI team reviews every complete submission blind, without names or schools.
- **Finals:** finalists present live to a panel of AI researchers. Five minutes, then questions.
- **Likely questions:** "How do you know it works?", "What would cause it to fail?", "Who would use this?"

**What stands out**

- Openness about limitations.
- A comparison with a simple baseline.
- A demo that runs reliably from start to finish.
- Reasons for each choice: why this problem, this model, these results.

**Common mistakes**

- Results too good to be true.
- A tutorial dataset with no new question.
- A write-up that lists tools but not what the team learned.
- Links judges cannot open.

## Awards

- An AI internship with a partner company for the Gold team.
- Scholarships to Veritas AI programs.
- Category awards: Best Beginner Build, Best AI for Good, Judges' Choice, and Merit Awards.
- A certificate for every participant who submits.

---

## How BRCA Spellcheck maps to the rubric

| Requirement | How the project meets it |
| --- | --- |
| Specific user | A genetic counselor with a patient's uncertain BRCA1 result |
| Input and output | A mutation in; a verdict, a score, the lab result and a track-record line out |
| Baseline | "Always harmless" and a 1992 substitution table, scored beside the model |
| Unseen test data | A whole region of the protein held out, tested once |
| Missed cases and false alarms | Reported separately for every method |
| Where it fails | Ten named wrong cases, one explained, with the cause labelled as a guess |
| Groups check | By protein region and by kind of swap. The data are mutations, not people, so age and sex cannot be checked |
| Public data | Lab measurements on cells; no patient records |
| Disclaimer | In the write-up, the notebook and the demo |
| Creativity | A protein language model graded against CRISPR lab truth |
| Clarity | One takeaway sentence, one table, one heat map, one typed mutation |
| Honest results | The false-alarm cost and the result that labels did not help are both stated |
